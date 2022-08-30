
const delaunay_chackbox = document.getElementById("delaunay");
const refresh_button = document.getElementById("refresh");
const submit_button = document.getElementById("submit-values");

delaunay_chackbox.addEventListener('change', (e) => {

    if (show_triang == 0) show_triang = 1
    else show_triang = 0;
    
    if (e.target.checked) {
        link.attr("stroke", (d) => d.etype[0] === "border" ? "black" : d.etype[0] === "rigid" ? "blue" : "red");
    } 
    else {
        link.attr('stroke', d => d.etype[0] === "border" ? "black" : null);
    }
});

refresh_button.addEventListener("click", () => {
    link = squishedCircle();
});

submit_button.addEventListener("click", () => {
    let values = document.getElementById("values").value;
    values = values.split(",");
    values = values.map((value, i) => [i, parseInt(value)]);
    link = squishedCircle(values);
});

