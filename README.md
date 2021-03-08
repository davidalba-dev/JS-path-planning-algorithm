# JS-path-planning-algorithm
# How to test
1) install http-server node module globally.
2) run the folowing command in terminal without quote in the project folder.
    "http-server"
3) run command "node dron-path.js" without quote.
4) get to "localhost:8080/show.html" in chrome.

# Input format
First line contains n, index, angle, interLineAge, direction

n means the number of polygon nodes.
index is the index of entry point in polygon nodes list - [1, n]
interLineAge is the distance between parallel lines
direction is 1 or -1, where 1 means go right otherwise, goes left.

Next n lines contain ordinates of polygon nodes on x-axis, y-axis in each line.
