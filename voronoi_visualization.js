var c = document.getElementById("myCanvas");
c.height = 500;
c.width = 500;
var ctx = c.getContext("2d");

function drawCell(text) {
    points = text.split(" ");
    ctx.beginPath();
    ctx.moveTo(10*Number(points[0]), 10*Number(points[1]));
    for (var i = 2; i < points.length-1; i = i + 2) {
        ctx.lineTo(10*Number(points[i]), 10*Number(points[i+1]));
    }
    ctx.closePath();
    ctx.stroke();
}

function drawCells(text) {
    var lines = text.split("\n");

    for (var i = 0; i < lines.length; i++)
    {
        drawCell(lines[i]);
    }
}

// https://stackoverflow.com/a/59804354/12724988
async function loadFile(file) {
    let text = await file.text();
    drawCells(text);
}