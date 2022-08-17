
let [w, h, s] = Array(3).fill(450);

let svg = d3.select('body')
            .append('svg')
            .attr('width', w)
            .attr('height', h);

const xScale = d3.scaleLinear()
                 .domain([0, 1])
                 .range([0, w]);
            
const yScale = d3.scaleLinear()
                 .domain([0, 1])
                 .range([h, 0])


function poissonSampling (n, width, height) {
    let d = width / Math.sqrt((n * height) / width);
    let pds = new PoissonDiskSampling({
        shape: [width, height],
        minDistance: d * 0.8,
        maxDistance: d * 1.6,
        tries: 15
    });
    return pds.fill();
}




// let data = [["A", 10], ["B", 100], ["C", 50], ["D", 5]];
// let data = [["A",10], ["B", 100], ["C", 50], ["D",5], ["E",1]];
// let data = [["A", 10], ["B", 100], ["C", 50], ["D", 5], ['E', 300], ['F', 50],  ['G', 400]];

let max_value = 500;
let max_circles = 20;
let n_data = Math.floor(Math.random() * max_circles) + 1;
let data = [...Array(n_data)].map((i) => [i, Math.floor(Math.random() * max_value)]);

let hier = d3.hierarchy(
    {name: "root", children: data.map( ([name, value]) => ({name, value}) )}
).sum(d => d.value);

let packed = d3.pack().padding(0.01)(hier);


// Polygonize circles

const edgeSz = 0.04;
let vertices = [];
let constraints = []; // Ex: [i, j] => v[i] and v[j] form a constraint
let circles = [];
let vtype = []; // Ex: vtype[i] == ['border', 0] => v[i] is the border of circle 0
let n_borders = [];

for (let d of [packed, ...packed.children]) {
    let {x, y, r} = d;
    const n = Math.max(12, Math.round(2*Math.PI * r / edgeSz)); // number of points in the border of each circle
    n_borders.push(n);
    let prev = vertices.length + n - 1;
    let icircle = circles.length;
    circles.push({x: xScale(x), y: yScale(y), r: xScale(r), ivtx: vertices.length, nvtx: n})

    for (let i = 0; i < n; i++) {
        const ang = 2*Math.PI * i / n;
        let px = xScale(r * Math.cos(ang) + x);
        let py = yScale(r * Math.sin(ang) + y);
        constraints.push([prev, vertices.length]);
        prev = vertices.length;
        vertices.push([px, py]);
        vtype.push(["border", icircle]);
    }
}


// Add a few random vertices

let {x, y, r} = circles[0]; // The outer circle
let empty_space = false;

for (let p of poissonSampling(800, 1, 1)) {
    p =[xScale(p[0]), yScale(p[1])];
    if (Math.hypot(p[0] - x, p[1] - y) > r*0.96) continue; // Outside outer circle
    let icircle = 0;  
    for (let i = 1; i < circles.length; i++) {
        let {x, y, r} = circles[i];
        if (Math.hypot(p[0] - x, p[1] - y) < r) {
            icircle = i;
            break;
        }
    }

    if (empty_space) {
        if (icircle > 0) {
            vertices.push(p);
            vtype.push(["internal", icircle]);
        }
    }
    else {
        vertices.push(p);
        vtype.push(["internal", icircle]);
    }
}


// compute the inner circles desired size

let totalArea = 0;
for (let c of circles) {
    c.area = Math.PI * c.r**2;
    totalArea += c.area;
}
let hullArea = circles[0].area;
let innerArea = totalArea - hullArea;
let areaEnlargement = hullArea / innerArea;
let newInnerArea = 0;
for (let c of circles.slice(1)) {
    c.newArea = c.area * areaEnlargement;
    c.newR = Math.sqrt(c.newArea / Math.PI);
    newInnerArea += c.newArea
}


// Add additional points and create a constrained Delaunay triangulation

// Creates a constrained triangulation, and populates the mesh with fields edges and etypes

const del = d3.Delaunay.from(vertices);
const {points, triangles, halfedges} = del;

const con = new Constrainautor(del._delaunator);
const conEdge = [];

for (const [v1, v2] of constraints) {
    con.constrainOne(v1, v2);
    conEdge[v1] = v2;
}
con.delaunify();


let edges = [];
let etype = [];

for (let i = 0, n = halfedges.length; i < n; ++i) {
    let j = halfedges[i];
    if (j < i) continue;
    const vi = triangles[i];
    const vj = triangles[j];
    edges.push([vi, vj]);

    if (vtype[vi][0] === "border" && vtype[vj][0] === "border") {
        if (vtype[vi][1] == vtype[vj][1]) {
            etype.push( [conEdge[vi] == vj || conEdge[vj] == vi ? "border" : "rigid", vtype[vi][1]] );
        }
        else etype.push(["shrink", 0]);
    }
    else {
        let ivertex = vtype[vi][0] == "internal" ? vtype[vi][1] : vtype[vj][1];
        if (ivertex == 0) etype.push(["shrink", 0])
        else etype.push(["rigid", ivertex])
    }
}


const nhull = circles[0].nvtx;
let node_data = vertices.map(([x, y], i) => {
    let obj = { x, y, vtype : vtype[i] };
    if (i < nhull) {
        [obj.fx, obj.fy] = [obj.x, obj.y];
    }
    return obj;
})

let links = edges.map(([isrc, idst], iedge) => {

    let obj = { source: node_data[isrc], target: node_data[idst], etype: etype[iedge] };
    obj.esize = Math.hypot(obj.source.x - obj.target.x, obj.source.y - obj.target.y);

    if (obj.etype[0] === "rigid") {
        let icircle = obj.etype[1];
        obj.dist = obj.esize * circles[icircle].newR / circles[icircle].r;
    }
    else {
        obj.dist = 0;
    }
    return obj
});

let links_copy = links.map(d => Object.create(d));






function paint_area(node_data, c_node_data, n_borders, c_i) {
    let path = d3.path();
    let init = d3.sum(n_borders.slice(0, c_i + 1));
    let n = d3.sum(n_borders.slice(0, c_i + 2));
    for (let i=init; i<n; i++) {
        path.moveTo(node_data[i].x, node_data[i].y)
        if (i != n-1) path.lineTo(node_data[i+1].x, node_data[i+1].y)
        else path.lineTo(node_data[init].x, node_data[init].y)
        path.lineTo(c_node_data[c_i].x, c_node_data[c_i].y)
    }
    return path;
}

function center_polygon(node_data, n_borders, c_i) {
    let init = d3.sum(n_borders.slice(0, c_i));
    let n = d3.sum(n_borders.slice(0, c_i+1));
    let n_border = n_borders[c_i];
    let c_coord = { x: 0, y: 0};

    for (let i=init; i<n; i++) {
        c_coord.x += node_data[i].x / n_border;
        c_coord.y += node_data[i].y / n_border;
    }
    return c_coord;
}





let c_node_data = [...Array(circles.length-1)].map((d, i) => center_polygon(node_data, n_borders, i+1));

let c_j = 3;
let path = paint_area(node_data, c_node_data, n_borders, c_j);

// let colors = ['gray', 'blue', 'green', 'yellow', 'red'];
let colors = ['#588c7e', '#f2e394', '#f2ae72', '#d96459'];
let paths = [...Array(n_borders.length-1)].map((d, i) => paint_area(node_data, c_node_data, n_borders, i))





svg.selectAll("line")
   .data(constraints.slice(0, circles[0].nvtx))
   .enter()
   .append("line")
   .attr('x1', d => vertices[d[0]][0])
   .attr('y1', d => vertices[d[0]][1])
   .attr('x2', d => vertices[d[1]][0])
   .attr('y2', d => vertices[d[1]][1])
   .attr("stroke", "black");

let areas = svg.selectAll('path')
               .data(paths)
               .enter()
               .append('path')
               .attr("d", (d, i) => paths[i])
               .attr('opacity', 0.5)
               .attr("fill", (d, i) => colors[i%colors.length]);

const link = svg.selectAll(".link")
                .data(links)
                .enter()
                .append("line")
                .classed("link", true)
                .classed("border", (d) => d.etype[0] === "border")
                .attr('x1', d => d.source.x)
                .attr('y1', d => d.source.y)
                .attr('x2', d => d.target.x)
                .attr('y2', d => d.target.y)
                // .attr("stroke", (d) => d.etype[0] === "border" ? "black" : d.etype[0] === "rigid" ? "blue" : "red")
                .attr('stroke', d => d.etype[0] === "border" ? "black" : null)
                // .attr('stroke', d => "black")

const node = svg.selectAll(".node")
                .data(node_data)
                .enter()
                .append("circle")
                .classed("node", true)
                .classed("fixed", (d) => d.fx !== undefined)
                .attr('cx', d => d.x)
                .attr('cy', d => d.y)
                .attr("r", 4)
                .attr("fill", "black")
                .attr("opacity", 0)

const c_node = svg.selectAll(".c_node")
                .data(c_node_data)
                .enter()
                .append("circle")
                .attr('cx', d => d.x)
                .attr('cy', d => d.y)
                // .attr("r", (d, i) => 5)
                // .attr("opacity", 0)

const text = svg.selectAll("text")
                .data(c_node_data)
                .enter()
                .append("text")
                .text((d, i) => data[i][1])
                .attr("x", (d, i) => c_node_data[i].x)
                .attr("y", (d, i) => c_node_data[i].y + circles[i+1].r * 0.25)
                .attr("font-family", "sans-serif")
                .attr("font-size", (d, i) => circles[i+1].r * 0.7)
                .style("text-anchor", "middle")



function dist(u, v=[0, 0]) {
    return Math.sqrt( Math.pow((u[0]-v[0]), 2) + Math.pow((u[1]-v[1]), 2) );
}

function dot_prod(u, v) {
    return u[0] * v[0] + u[1] * v[1];
}

function triang_heights(points) {
    let heights = []
    for (let i=0; i<3; i++) {
        let [x_0, y_0] = points[i];
        let [x_1, y_1] = points[(i+1)%3];
        let [x_2, y_2] = points[(i+2)%3];
        // let base = dist([x_1, y_1], [x_2, y_2]);
        let v = [x_2 - x_1, y_2 - y_1];
        let u = [x_0 - x_1, y_0 - y_1];
        let proj = dot_prod(u, v) / (dist(v)**2);
        let u_proj = [u[0] - proj*v[0], u[1] - proj*v[1]];
        heights.push(dist(u_proj));
    }
    return heights;
}

function perpendiculars(points) {
    let perps = [];
    for (let i=0; i<3; i++) {
        let [x_1, y_1] = points[(i+1)%3];
        let [x_2, y_2] = points[(i+2)%3];
        let v = [x_2 - x_1, y_2 - y_1];
        let mod_v = dist(v);
        v = [v[0]/mod_v, v[1]/mod_v];
        perps.push([v[1], -v[0]]);
    }
    return perps;
}

const check_border = (vtype) => vtype == 'border';

let n_triang = triangles.length/3 -1;
let circles_t = [...Array(circles.length)].map(() => 0);

for (let i=0; i<n_triang; i++) {
    let t = [...Array(3)].map((d, j) => node_data[triangles[3*i + j]].vtype);

    if (!t.some(check_border)) 
        circles_t[t[0][1]] += 1;
    else {
        for (let j=0; j<3; j++) {
            if (t[j][0] != 'border') circles_t[t[j][1]] += 1;
        }
    }
}

let x_0, y_0, x_1, y_1, x_2, y_2, v_x, v_y, p_v, d, b_x, b_y, dot, bar_d;
let cte = 0.5;
let d_dot = 0.1;
let delta_d = 5;
let scale = 50;
let apply_force = true;
// apply_force = false;
let sv;
let c_i, A, delta_A, delta_A_max;

let sim = d3.forceSimulation(node_data)
            .force("link", d3.forceLink(links)
                             .distance(d => d.dist)
                             .strength((d) => (d.dist == 0 ? 1 : 0.03)))
            .force('triangle-force', () => {

                if (apply_force) {
                    for (let i=0; i < n_triang; i++) {

                        let types = [...Array(3)].map((d, j) => node_data[triangles[3*i + j]].vtype[0]);
                        let vtypes = [...Array(3)].map((d, j) => node_data[triangles[3*i + j]].vtype);

                        for (let j=0; j<3; j++) {
                            if (vtypes[j][0] != 'border') c_i = vtypes[j][1]
                        }
                        
                        if (c_i != 0 && !types.every(check_border)) {
                            A = circles[c_i].newArea / circles_t[c_i];

                            if (!types.some(check_border)) {
                                s = 0.4;
                                sv = s;
                                delta_A = A * 0.5;
                                delta_A_max = 0;
                            }
                            else {
                                s = 1.1;
                                sv = 0.5;
                                delta_A = A * 0.9;
                                delta_A_max = A * 0.3
                            }

                            [x_0, y_0] = [node_data[triangles[3*i]].x, node_data[triangles[3*i]].y];
                            [x_1, y_1] = [node_data[triangles[3*i+1]].x, node_data[triangles[3*i+1]].y];
                            [x_2, y_2] = [node_data[triangles[3*i+2]].x, node_data[triangles[3*i+2]].y];
        
                            d = [
                                (x_0 - x_1) * (y_2 - y_1) - (y_0 - y_1) * (x_2 - x_1),
                                (x_1 - x_2) * (y_0 - y_2) - (y_1 - y_2) * (x_0 - x_2),
                                (x_2 - x_0) * (y_1 - y_0) - (y_2 - y_0) * (x_1 - x_0)
                            ];
        
                            let base = dist([x_1, y_1], [x_2, y_2]);
                            let heights = triang_heights([[x_0, y_0], [x_1, y_1], [x_2, y_2]]);
                            let Area = base * heights[0] / 2;
        
                            
        
                            let min_i = heights.indexOf(Math.min(...heights));
                            let max_i = heights.indexOf(Math.max(...heights));
                            let perps = perpendiculars([[x_0, y_0], [x_1, y_1], [x_2, y_2]]);
        
                            if (d[min_i] > 0) {
                                if (Area >= A + delta_A + delta_A_max) {
                                    node_data[triangles[3*i+max_i]].x -= s * perps[max_i][0];
                                    node_data[triangles[3*i+max_i]].y -= s * perps[max_i][1];
                                    node_data[triangles[3*i+max_i]].vx -= sv * perps[max_i][0];
                                    node_data[triangles[3*i+max_i]].vy -= sv * perps[max_i][1];
                                }
                                if (Area < A - delta_A) {
                                    node_data[triangles[3*i+min_i]].x += s * perps[min_i][0];
                                    node_data[triangles[3*i+min_i]].y += s * perps[min_i][1];
                                    node_data[triangles[3*i+min_i]].vx += sv * perps[min_i][0];
                                    node_data[triangles[3*i+min_i]].vy += sv * perps[min_i][1];
                                }
                            }
                            
                            if (d[min_i] < 0) {
                                node_data[triangles[3*i+min_i]].x += s * perps[min_i][0];
                                node_data[triangles[3*i+min_i]].y += s * perps[min_i][1];
                                node_data[triangles[3*i+min_i]].vx += sv * perps[min_i][0];
                                node_data[triangles[3*i+min_i]].vy += sv * perps[min_i][1];
                            }
                        }
                    }
                }
            })



function dragged(event, d) {
    d.x = event.x;
    d.y = event.y;
    sim.alpha(0.3).restart();
}

const drag = d3.drag()
               .on('drag', dragged)


sim.on('tick', () => {
    node.attr('cx', d => d.x)
        .attr('cy', d => d.y)
        .call(drag);

    paths = paths.map((d, i) => paint_area(node_data, c_node_data, n_borders, i))
    areas.attr("d", (d, i) => paths[i])

    link.attr('x1', d => d.source.x)
    .attr('y1', d => d.source.y)
    .attr('x2', d => d.target.x)
    .attr('y2', d => d.target.y)

    c_node_data = c_node_data.map((d, i) => center_polygon(node_data, n_borders, i+1));

    c_node.attr('cx', (d, i) => c_node_data[i].x)
            .attr('cy', (d, i) => c_node_data[i].y);

    text.attr("x", (d, i) => c_node_data[i].x)
        .attr("y", (d, i) => c_node_data[i].y + circles[i+1].r * 0.25)
        .style("text-anchor", "middle");
});