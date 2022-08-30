
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

function dragged(event, d) {
    d.x = event.x;
    d.y = event.y;
}


function forceTriangle(triangles = []) {
    let area = (t, a0, aCurrent) => a0; // Default function that maps triplet to triangle area
    let strength = (t) => 1; // Default function that maps triplet to force strength
    let iterations = 1; // Default number of iterations
    let nodes; // Array of nodes
    let areas; // Array of initial triangle areas
    const makeTriangle = (
        triplet // Creates a Triangle object from a node triplet
    ) =>
        new Triangle(
        ...triplet.map((node) => [node.x + node.vx, node.y + node.vy])
        );

    function force(alpha) {

        for (let k = 0; k < iterations; ++k) {
            for (let itriangle = 0; itriangle < triangles.length; itriangle++) {

                let triplet = triangles[itriangle];
                let a0 = areas[itriangle];
                let t = makeTriangle(triplet);
                let q = t.p.map(([x, y]) => [x, y]); // Clone vertices
                let newArea = area(triplet, a0, t.area());
                t.setArea(newArea);
                let s = strength(triplet);
                for (let i = 0; i < 3; i++) {
                    let node = triplet[i];
                    node.vx += (t.p[i][0] - q[i][0]) * alpha * s;
                    node.vy += (t.p[i][1] - q[i][1]) * alpha * s;
                }
            }
        }
    }

    function initialize(_nodes, _random) {
        nodes = _nodes;
        areas = [];
        for (let triplet of triangles) {
            areas.push(makeTriangle(triplet).area());
        }
    }

    force.initialize = initialize;

    force.triangles = function (_) {
        return arguments.length ? ((triangles = +_), force) : triangles;
    };

    force.iterations = function (_) {
        return arguments.length ? ((iterations = +_), force) : iterations;
    };

    force.strength = function (_) {
        return arguments.length
        ? ((strength = typeof _ === "function" ? _ : () => _), force)
        : strength;
    };

    force.area = function (_) {
        return arguments.length
        ? ((area = typeof _ === "function" ? _ : () => _), force)
        : area;
    };

    return force;
}


class Triangle {
    constructor(a, b, c) {
        this.p = [a, b, c];
    }

    area() {
        let total = 0;
        let [px, py] = this.p[2];
        for (let [x, y] of this.p) {
            total += (px - x) * (y + py);
            [px, py] = [x, y];
        }
        return total / 2;
    }

    sideLength(i) {
        const j = (i + 1) % 3;
        return Math.hypot(this.p[i][0] - this.p[j][0], this.p[i][1] - this.p[j][1]);
    }

    sideHeight(i) {
        const j = (i + 1) % 3,
        k = (i + 2) % 3;
        const u = [this.p[j][0] - this.p[i][0], this.p[j][1] - this.p[i][1]];
        const v = [this.p[k][0] - this.p[j][0], this.p[k][1] - this.p[j][1]];
        const ulen = Math.hypot(...u) || 0.0001;
        const uhat = [u[0] / ulen, u[1] / ulen];
        const uhatDotV = uhat[0] * v[0] + uhat[1] * v[1];
        return [v[0] - uhatDotV * uhat[0], v[1] - uhatDotV * uhat[1]];
    }

    barycenter() {
        const sum = this.p.reduce((a, b) => [a[0] + b[0], a[1] + b[1]]);
        return [sum[0] / 3, sum[1] / 3];
    }

    setArea(a) {
        let factor = a / this.area();
        if (Math.abs(factor) < 1) return this.setArea1(factor);
        return this.setArea2(factor);
    }

    setArea1(factor) {
        // Use scaling towards barycenter
        //let formerArea = this.area();
        let f = Math.sqrt(Math.abs(factor));
        let bary = this.barycenter();
        let v = [];
        this.p.forEach((q) => {
            let u = [(q[0] - bary[0]) * f, (q[1] - bary[1]) * f];
            v.push(u);
            q[0] = bary[0] + u[0];
            q[1] = bary[1] + u[1];
        });
        if (factor < 0) {
            // Reflect wrt barycenter using the smallest height as an axis
            let { i, h, hgt } = [0, 1, 2]
                .map((i) => {
                    let hgt = this.sideHeight(i);
                    let h = Math.hypot(...hgt);
                    return { i, h, hgt };
                })
                .reduce((a, b) => (a.h < b.h ? a : b));

            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            let n = [hgt[0] / h, hgt[1] / h];
            let a = v[k][0] * n[0] + v[k][1] * n[1];
            let b = h - a;
            this.p[i] = [this.p[i][0] + n[0] * b * 2, this.p[i][1] + n[1] * b * 2];
            this.p[j] = [this.p[j][0] + n[0] * b * 2, this.p[j][1] + n[1] * b * 2];
            this.p[k] = [this.p[k][0] - n[0] * a * 2, this.p[k][1] - n[1] * a * 2];
            // mutable debug = { hgt, ratio: this.area() / formerArea };
        }
    }

    setArea2(factor) {
        const h = [0, 1, 2].map((i) => this.sideHeight(i));
        const hlen = h.map(([x, y]) => Math.hypot(x, y));
        let i =
        Math.abs(factor) < 1
            ? // Shrink
                hlen[0] > hlen[1]
                ? hlen[0] > hlen[2]
                    ? 0
                    : 2
                : hlen[1] > hlen[2]
                ? 1
                : 2
            : // Expand
            hlen[0] < hlen[1]
            ? hlen[0] < hlen[2]
                ? 0
                : 2
            : hlen[1] < hlen[2]
            ? 1
            : 2;
        let u = h[i],
            v = [0, 0];
        if (factor < 0) {
            factor = -factor;
            v = [...u];
            u = [-u[0], -u[1]];
        }
        u = [v[0] + (-u[0] * (factor - 1)) / 3, v[1] + (-u[1] * (factor - 1)) / 3];
        const j = (i + 1) % 3,
              k = (i + 2) % 3;
        this.p[i][0] += u[0];
        this.p[i][1] += u[1];
        this.p[j][0] += u[0];
        this.p[j][1] += u[1];
        this.p[k][0] -= u[0] * 2;
        this.p[k][1] -= u[1] * 2;
    }
}
