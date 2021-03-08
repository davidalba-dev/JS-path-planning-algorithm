'use strict';

const fs = require('fs');
const ws = fs.createWriteStream("out.txt");
const ep = 1.0e-7;

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

function crossProduct(A, B) { return A.scala(B); }

/// Get Convex Hull -- Graham's Algorithm ///
function getConvexHull(poly) {
    poly.sort(function(A, B) {
        if (A.x != B.x) return A.x - B.x;
        return A.y - B.y;
    });

    /// erase same points and make polygon list unique ///
    let filteredPoly = [];
    for (let i = 0; i < poly.length; i ++) {
        if (i == 0) {
            filteredPoly.push(poly[0]);
        } else {
            if (poly[i].x - poly[i - 1].x > ep || poly[i].y - poly[i - 1].y > ep) {
                filteredPoly.push(poly[i]);
            }
        }
    }

    poly = filteredPoly;
    let n = poly.length;
    let m = 0;
    let ch = [];

    for (let i = 0; i < n; i ++) {
        while (m > 1 && crossProduct(ch[m - 1].minus(ch[m - 2]), poly[i].minus(ch[m - 2])) <= 0.0) m --;
        ch[m ++] = poly[i];
    }

    let k = m;
    for (let i = n - 2; i >= 0; i --) {
        while (m > k && crossProduct( ch[m - 1].minus(ch[m - 2]), poly[i].minus(ch[m - 2]) ) <= 0.0) m --;
        ch[m ++] = poly[i];
    }
    if (n > 1) m --;
    ch.length = m;
    return ch;
}

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

function insideConvexPoly(poly, q) {
    let cnt = 0;
    let n = poly.length;
    for (let i = 0; i < n; i ++) {
        let j = (i + 1) % n;
        if (inside(q, poly[i], poly[j])) return 2;
        if (poly[i].minus(q).scala(poly[j].minus(q)) > ep) {
            cnt ++;
        }
    }
    if (cnt == n || !cnt) return 1;
    return 0;
}

function getIntersectionSegmentsInsideConvex(poly, entryPoint, angle, interLineAge, direction) {
    let rad = angle * Math.PI / 180;
    let dx, dy;
    dx = 1.0 * Math.cos(rad);
    dy = Math.sin(rad);
    // if (Math.abs(rad) - Math.PI / 2 < ep) dy = 0;
    // else dy = Math.tan(rad);

    let n = poly.length;

    for (let step = 0; ; step ++) {
        let A, B;

        if (direction == 1) { /// right direction
            if (Math.sin(rad) < ep) {
                A = entryPoint.plus(new pnt(0.0, -interLineAge).multi(parseFloat(step)));
                B = A.plus(new pnt(dx, dy));
            } else {
                A = entryPoint.plus(new pnt(interLineAge / Math.sin(rad), 0.0).multi(parseFloat(step)));
                B = A.plus(new pnt(dx, dy));
            }
        } else if (direction == -1) { /// left direction
            if (Math.sin(rad) < ep) {
                A = entryPoint.minus(new pnt(0.0, -interLineAge).multi(parseFloat(step)));
                B = A.plus(new pnt(dx, dy));
            } else {
                A = entryPoint.minus(new pnt(interLineAge / Math.sin(rad), 0.0).multi(parseFloat(step)));
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
            if (A.x != B.x) return A.x - B.x;
            return A.y - B.y;
        });

        let filteredIntersections = [];
        for (let i = 0; i < intersections.length; i ++) {
            if (i == 0) {
                filteredIntersections.push(intersections[0]);
            } else {
                if (intersections[i].x - intersections[i - 1].x > ep || intersections[i].y - intersections[i - 1].y > ep) {
                    filteredIntersections.push(intersections[i]);
                }
            }
        }
        intersections = filteredIntersections;

        for (let i = 0; i < intersections.length - 1; i ++) {
            let overlapping = false;
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

            if (overlapping == true) continue;

            /// choose segments inside polygon ---- for double-check ///
            let middlePoint = intersections[i].plus(intersections[i + 1]).divide(2);
            let insideConvex = insideConvexPoly(poly, middlePoint);
            if (insideConvex) {
                if (intersections[i].y < intersections[i + 1].y) {
                    segments.push(new segment(intersections[i], intersections[i + 1])); /// seg.A : lower y ordinate
                } else {
                    segments.push(new segment(intersections[i + 1], intersections[i])); /// seg.B : upper y ordinate
                }
            }
        }
    }
}

function getIntersectionSegments(poly, entryPoint, angle, interLineAge, direction) {
    poly = getConvexHull(poly); /// anti-clock wised convex-hull polygon

    getIntersectionSegmentsInsideConvex(poly, entryPoint, angle, interLineAge, direction);

    /// setting drone-path ///
    let path = [];
    path.push(entryPoint);
    let id;

    if (segments.length) {
        if (entryPoint.distance(segments[0].A) + ep < entryPoint.distance(segments[0].B)) {
            path.push(segments[0].A);
            path.push(segments[0].B);
            id = 1;
        } else {
            path.push(segments[0].B);
            path.push(segments[0].A);
            id = 0;
        }
    
        if (segments.length > 1) for (let i = 1; i < segments.length; i ++) {
            if (id) {
                path.push(segments[i].B);
                path.push(segments[i].A);
            } else {
                path.push(segments[i].A);
                path.push(segments[i].B);
            }
            id = 1 - id;
        }
    }

    /// output ///
    ws.write(poly.length + '\n');
    for (let i = 0; i < poly.length; i ++) {
        ws.write(poly[i].x + " " + poly[i].y + "\n");
    }

    ws.write(path.length + '\n');
    for (let i = 0; i < path.length; i ++) {
        ws.write(path[i].x + " " + path[i].y + "\n");
    }
    ws.end();
}

function main() {
    fs.readFile('in.txt', 'utf-8', (err, data) => {
        if (err) throw err; 
      
        data = data.replace(/\s+$/g, '').split('\r');
        let spec = data[0].split(' ');
        const n = parseInt(spec[0], 10);
        const entry = parseInt(spec[1], 10);
        const angle = parseFloat(spec[2]);
        const interLineAge = parseFloat(spec[3]);
        const direction = parseInt(spec[4]);

        let poly = [];

        for (let i = 1; i <= n; i ++) {
            const latLong = data[i].split(' ').map(temp => parseFloat(temp));
            poly.push(new pnt(latLong[0], latLong[1]));
        }
        getIntersectionSegments(poly, poly[entry - 1], angle, interLineAge, direction);
    });
}

main();