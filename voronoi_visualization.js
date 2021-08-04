var c = document.getElementById("myCanvas");
var ctx = c.getContext("2d");

function drawCell(text) {
    points = text.split(" ");

    ctx.fillStyle = "rgb(" + points[0] + "," + points[1] + "," + points[2] + ")";

    ctx.beginPath();
    ctx.moveTo(Number(points[3]), Number(points[4]));

    for (var i = 5; i < points.length-1; i = i + 2) {
        ctx.lineTo(Number(points[i]), Number(points[i+1]));
    }
    
    ctx.closePath();
    ctx.fill();
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

async function updateDimensions(file) {
    let text = await file.text();
    var firstLine = text.split("\n")[0];
    c.width = Number(firstLine.split(" ")[0]);
    c.height = Number(firstLine.split(" ")[1]);
}