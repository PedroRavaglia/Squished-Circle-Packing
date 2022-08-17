
const delaunay_chackbox = document.getElementById("delaunay");

delaunay_chackbox.addEventListener('change', (e) => {
    if (e.target.checked) {
        link.attr("stroke", (d) => d.etype[0] === "border" ? "black" : d.etype[0] === "rigid" ? "blue" : "red");
    } 
    else {
        link.attr('stroke', d => d.etype[0] === "border" ? "black" : null);
    }
});