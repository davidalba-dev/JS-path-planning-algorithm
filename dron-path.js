'use strict';

const fs = require('fs');
const ws = fs.createWriteStream("out.txt");

const entryPointTitleList = ["topLeft", "topRight", "bottomLeft", "bottomRight"];
const offsets = [1, 1, -1, -1];

const LAT_LONG_MODE = true;

let segments = [];

function pnt(x, y) {
    this.x = x;
    this.y = y;
}

function segment(A, B) {
    this.A = A;
    this.B = B;
}

pnt.prototype.minus = function(A) {
    return new pnt(this.x - A.x, this.y - A.y);
};

pnt.prototype.plus = function(A) {
    return new pnt(this.x + A.x, this.y + A.y);
};

pnt.prototype.multi = function(k) {
    return new pnt(this.x * parseFloat(k), this.y * parseFloat(k));
};

pnt.prototype.divide = function(k) {
/// k -> 0.0 return null ///
    return new pnt(this.x / parseFloat(k), this.y / parseFloat(k));
};

pnt.prototype.distance = function (A) {
    let P = this.minus(A);
    return Math.sqrt(P.vector(P));
}

pnt.prototype.scala = function(A) {
    return this.x * A.y - this.y * A.x;
};

pnt.prototype.vector = function(A) {
    return this.x * A.x + this.y * A.y;
};

function convertMeterToLatLong(distanceM){
    return (0.0001 / 1.11) * distanceM
}

const ep = LAT_LONG_MODE ? convertMeterToLatLong(1.0e-9) : 1.0e-9;

function inside(A, B, C) { /// if pnt A is on the segment between pnt B and pnt C, return 1
    let a = A.minus(B).vector(C.minus(B));
    let b = C.minus(B).vector(C.minus(B));
    let c = A.minus(B).scala(C.minus(B));
    if (Math.abs(c) > ep) return 0;
    if (a > -ep && a < b + ep) return 1;
    return 0;
}

/// intersect of two lines : Line (A, B), Segment (C, D) return intersection point ///
function intersect(A, B, C, D) {
    let den = A.minus(B).scala(C.minus(D)); /// (a - b) ^ (c - d)

	if (Math.abs(den) < ep) {		//The two lines are parallel...
        if (Math.abs(C.minus(B).scala(C.minus(A))) < ep) { /// (c - b) ^ (c - a)
            return null; /// overlapping
        }
        return null; /// parallel
    }
    let x = A.minus(D).scala(C.minus(D)) / den; /// (a - d ^ c - d) / den
    return A.multi(1 - x).plus(B.multi(x)); /// a * (1 - x) + b * x   The two lines are intersect...
}

function getCentroidPolygon(poly) {
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    let n = poly.length;
    for (let i = 0; i < n; i ++) {
        minX = Math.min(minX, poly[i].x);
        maxX = Math.max(maxX, poly[i].x);
        minY = Math.min(minY, poly[i].y);
        maxY = Math.max(maxY, poly[i].y);
    }
    return new pnt(minX, minY).plus(new pnt(maxX, maxY)).divide(2);
}

function getIntersectionSegmentsInsideConvex(poly, centroid, angle, distanceBetweenLines, direction) {
    let rad = angle * Math.PI / 180;
    let dx, dy;
    let delta = LAT_LONG_MODE ? convertMeterToLatLong(1.0) : 1.0;
    dx = delta * Math.cos(rad);
    dy = delta * Math.sin(rad);

    let n = poly.length;

    for (let step = (direction == 1 ? 0 : 1); ; step ++) {
        let A, B;

        if (direction == 1) { /// right direction
            if (Math.sin(rad) < ep) {
                A = centroid.plus(new pnt(0.0, -distanceBetweenLines).multi(parseFloat(step)));
                B = A.plus(new pnt(dx, dy));
            } else {
                A = centroid.plus(new pnt(distanceBetweenLines / Math.sin(rad), 0.0).multi(parseFloat(step)));
                B = A.plus(new pnt(dx, dy));
            }
        } else if (direction == -1) { /// left direction
            if (Math.sin(rad) < ep) {
                A = centroid.minus(new pnt(0.0, -distanceBetweenLines).multi(parseFloat(step)));
                B = A.plus(new pnt(dx, dy));
            } else {
                A = centroid.minus(new pnt(distanceBetweenLines / Math.sin(rad), 0.0).multi(parseFloat(step)));
                B = A.plus(new pnt(dx, dy));
            }
        }

        let intersections = [];
        for (let i = 0; i < n; i ++) {
            let intersection = null;
            intersection = intersect(A, B, poly[i], poly[(i - 1 + n) % n]);
            if (intersection) { /// has intersection
                if (inside(intersection, poly[i], poly[(i - 1 + n) % n])) {
                    intersections.push(intersection); /// intersection is not null
                }
            }
        }
        
        if (intersections.length == 0) break;
        if (intersections.length == 1) continue;

        intersections.sort(function(A, B) {
            if (Math.abs(A.x - B.x) < ep) return A.x + ep - B.x;
            return A.y + ep - B.y;
        });

        let overlapping = false;
        for (let i = 0; i < intersections.length - 1; i ++) {
            /// check if each line is overlapping with polygon edges
            for (let j = 0; j < n; j ++) {
                let k = (j + 1) % n;
                let A = intersections[i];
                let B = intersections[i + 1];
                let C = poly[j];
                let D = poly[k];

                let den = A.minus(B).scala(C.minus(D)); /// (a - b) ^ (c - d)

                if (Math.abs(den) < ep) {		//The two lines are parallel...
                    if (Math.abs(C.minus(B).scala(C.minus(A))) < ep) { /// (c - b) ^ (c - a)
                        overlapping = true;
                    }
                }
            }
        }

        if (intersections.length == 2 && overlapping == true) continue;

        let m = intersections.length;
        if (intersections[0].y + ep < intersections[m - 1].y) {
            segments.push(new segment(intersections[0], intersections[m - 1])); /// seg.A : lower y ordinate
        } else if (intersections[0].y > ep + intersections[m - 1].y) {
            segments.push(new segment(intersections[m - 1], intersections[0])); /// seg.B : higher y ordinate
        } else {
            if (intersections[0].x + ep < intersections[m - 1].x) {
                segments.push(new segment(intersections[0], intersections[m - 1])); /// seg.A : lower x ordinate
            } else {
                segments.push(new segment(intersections[m - 1], intersections[0])); /// seg.B : higher x ordinate
            }
        }
    }
}

function getPath(poly, angle, distanceBetweenLines, entryPointTitle) {
    let centroid = getCentroidPolygon(poly);

    let direction;
    if (angle > 180) {
        angle -= 180;
        direction = 1; /// upper segments right side of transect 0
        getIntersectionSegmentsInsideConvex(poly, centroid, angle, distanceBetweenLines, direction);
        segments.reverse();

        direction = -1; /// upper segments left side of transect 0
        getIntersectionSegmentsInsideConvex(poly, centroid, angle, distanceBetweenLines, direction);

        for (let i = 0; i < segments.length; i ++) {
            segments[i] = new segment(segments[i].B, segments[i].A);
        }
    } else {
        direction = -1; /// upper segments left side of transect 0
        getIntersectionSegmentsInsideConvex(poly, centroid, angle, distanceBetweenLines, direction);
        segments.reverse();
    
        direction = 1;
        getIntersectionSegmentsInsideConvex(poly, centroid, angle, distanceBetweenLines, direction);
    }

    /// setting drone-path ///
    let n = segments.length;
    let entryID = entryPointTitleList.findIndex(title => title == entryPointTitle) 
    let startID = entryID > 1 ? n - 1 : 0;

    let path = [];
    let id;

    if (n) {
        if (entryPointTitle.indexOf('Left') > -1) id = 0;
        else id = 1;
        for (let i = startID, j = 0; j < n; j ++, i += offsets[entryID]) {
            if (!id) {
                path.push(segments[i].A);
                path.push(segments[i].B);
            } else {
                path.push(segments[i].B);
                path.push(segments[i].A);
            }
            id = 1 - id;
        }
    }

    /// For show on browser - need to be removed ///
    if (!LAT_LONG_MODE) {
        ws.write(poly.length + '\n');
        for (let i = 0; i < poly.length; i ++) {
            ws.write(poly[i].x + " " + poly[i].y + "\n");
        }
    } else {
        ws.write(poly.length + '\n');
        for (let i = 0; i < poly.length; i ++) {
            ws.write((poly[i].x * 100000 - 671500) + " " + (poly[i].y * 100000 - 4649000) + "\n");
        }
    }
    /// End -- For show on browser - need to be removed ///
    return path;
}

function main() {
    fs.readFile('in.txt', 'utf-8', (err, data) => {
        if (err) throw err; 
      
        data = data.replace(/\s+$/g, '').split('\n');
        let spec = data[0].split(' ');
        const n = parseInt(spec[0], 10);
        const angle = parseFloat(spec[1]);
        const distanceBetweenLines = parseFloat(spec[2]);
        const entryPointTitle = spec[3].split('\r')[0];

        let poly = [];

        for (let i = 1; i <= n; i ++) {
            const latLong = data[i].split(' ').map(temp => parseFloat(temp));
            poly.push(new pnt(latLong[0], latLong[1]));
        }

        /// Final Answer - path ///
        let path = getPath(poly, angle, LAT_LONG_MODE ? convertMeterToLatLong(distanceBetweenLines) : distanceBetweenLines, entryPointTitle);
    
        /// For show on browser - Need to be removed ///
        if (LAT_LONG_MODE) {
            let centroid = getCentroidPolygon(poly);
            ws.write((centroid.x * 100000 - 671500) + " " + (centroid.y * 100000 - 4649000) + "\n");

            ws.write(path.length + '\n');
            for (let i = 0; i < path.length; i ++) {
                ws.write((path[i].x * 100000 - 671500) + " " + (path[i].y * 100000 - 4649000) + "\n");
            }
        } else {
            let centroid = getCentroidPolygon(poly);
            ws.write(centroid.x + " " + centroid.y + "\n");

            ws.write(path.length + '\n');
            for (let i = 0; i < path.length; i ++) {
                ws.write(path[i].x + " " + path[i].y + "\n");
            }
        }
        /// End -- For show on browser - Need to be removed ///
        ws.end();
    });
}

main();