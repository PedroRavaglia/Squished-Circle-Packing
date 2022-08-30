
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




let max_value = 500;
let max_circles = 20;

let init = 0;
let show_triang = 0;

function squishedCircle(values=[]) {
    if (init == 1) {
        svg.remove();
        svg = d3.select('body')
                .append('svg')
                .attr('width', w)
                .attr('height', h);
    }
    init = 1;

    let n_data, data;
    if (values.length == 0) {
        n_data = Math.floor(Math.random() * max_circles) + 1;
        data = [...Array(n_data)].map((i) => [i, Math.floor(Math.random() * max_value)]);
    }
    else {
        data = values;
    }

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
    let empty_space = true;

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


    // Add area constraints to make sure triangles do not flip
    let triplets = [];
    for (let i = 0, n = triangles.length; i < n; i += 3) {
        const v = [i, i + 1, i + 2].map((j) => node_data[triangles[j]]);
        const icircle = v[0].vtype[1];
        if (icircle > 0 && icircle == v[1].vtype[1] && icircle == v[2].vtype[1]) {
            v.areaRatio = (circles[icircle].newR / circles[icircle].r) ** 2;
            triplets.push(v);
        } 
        else {
            v.areaRatio = 0.01;
            triplets.push(v);
        }
    }


    let c_node_data = [...Array(circles.length-1)].map((d, i) => center_polygon(node_data, n_borders, i+1));

    let c_j = 3;
    let path = paint_area(node_data, c_node_data, n_borders, c_j);

    let colors = ['#588c7e', '#f2e394', '#f2ae72', '#d96459'];
    let paths = [...Array(n_borders.length-1)].map((d, i) => paint_area(node_data, c_node_data, n_borders, i));


    let n_triang = triangles.length/3 - 1;
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

    let sim = d3.forceSimulation(node_data)
                .force("link", d3.forceLink(links)
                                .distance(d => d.dist)
                                .strength((d) => (d.dist == 0 ? 1 : 0.03)))
                .force('triangle-force', forceTriangle(triplets).area((t, a0) => a0 * t.areaRatio));


    const drag = d3.drag()
                .on('drag', dragged)


    sim.on('tick', () => {
        node.attr('cx', d => d.x)
            .attr('cy', d => d.y)
            .call(drag);

        paths = paths.map((d, i) => paint_area(node_data, c_node_data, n_borders, i));
        areas.attr("d", (d, i) => paths[i]);

        link.attr('x1', d => d.source.x)
        .attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x)
        .attr('y2', d => d.target.y);

        c_node_data = c_node_data.map((d, i) => center_polygon(node_data, n_borders, i+1));

        c_node.attr('cx', (d, i) => c_node_data[i].x)
                .attr('cy', (d, i) => c_node_data[i].y);

        text.attr("x", (d, i) => c_node_data[i].x)
            .attr("y", (d, i) => c_node_data[i].y + circles[i+1].r * 0.25)
            .style("text-anchor", "middle");
    });



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
                    .attr('stroke', d => {
                        if (show_triang)
                            return d.etype[0] === "border" ? "black" : d.etype[0] === "rigid" ? "blue" : "red";
                        else
                            return d.etype[0] === "border" ? "black" : null;
                    })

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

    return link;
}

let link = squishedCircle();