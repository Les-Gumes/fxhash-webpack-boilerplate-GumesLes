// import * as d3 from "d3-voronoi"
// import poly2tri from "poly2tri"


export function transposeMatrix(matrix) {
  // Simple transposition
  var result = new Array(matrix[0].length);
	for (var i = 0; i < matrix[0].length; ++i) {
		result[i] = new Array(matrix.length);
	}
	for (var i = 0; i < matrix.length; ++i) {
		for (var j = 0; j < matrix[0].length; ++j) {
			result[j][i] = matrix[i][j];
		}
	}
	return result;
}

export let R = (a=1) =>  fxrand()*a
export let Rz = (a=1) =>  a*(R(2)-1)
export let RP = () => {
  let p = new Point(R(2)-1, R(2)-1, R(2)-1)
  return divP(p, norm2dP(p))
}
let L = (x, y)=>(x*x+y*y)**0.5;
export let M = (t, a, b)=> (1 - t) * a + t * b;
export let Mp = (t, a, b)=> new Point(
  (1 - t) * a.x + t * b.x,
  (1 - t) * a.y + t * b.y)

export function norm2d(vector) {
  return Math.sqrt(vector[0]**2 + vector[1]**2)
}

export function norm3d(vector) {
  return (vector[0]**2 + vector[1]**2 + vector[2]**2)**(1/2)
}

export function hslToHex(h, s, l) {
  l /= 100;
  const a = s * Math.min(l, 1 - l) / 100;
  const f = n => {
    const k = (n + h / 30) % 12;
    const color = l - a * Math.max(Math.min(k - 3, 9 - k, 1), -1);
    return Math.round(255 * color).toString(16).padStart(2, '0');   // convert to Hex and prefix "0" if needed
  };
  return `#${f(0)}${f(8)}${f(4)}`;
}

export function multiplyMatrices(m1, m2) {
    // Naive multiplication
    var result = [];
    for (var i = 0; i < m1.length; i++) {
        result[i] = [];
        for (var j = 0; j < m2[0].length; j++) {
            var sum = 0;
            for (var k = 0; k < m1[0].length; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

function CatmullRom(tho) {
  // Creating a CatmullRom Matrix accoring to
  // https://lucidar.me/en/mathematics/catmull-rom-splines/
  let Mat = [
    [0, 1, 0, 0],
    [-tho, 0, tho, 0],
    [2*tho, tho-3, 3-2*tho, -tho],
    [-tho, 2 - tho, tho - 2, tho]
  ]
  return Mat
}

export function rotationMatrix(theta) {
  return [
    [Math.cos(theta), -Math.sin(theta)],
    [Math.sin(theta), Math.cos(theta)]
  ]
}

export function rotate2center(p, theta) {
    p.x = p.x * Math.cos(theta) - p.y * Math.sin(theta)
    p.y = p.x * Math.sin(theta) + p.y * Math.cos(theta)
    return p
}

export function rotate2d(point, theta, center=[0, 0]) {
  // Transpose according to center
  point[0] = point[0] - center[0]
  point[1] = point[1] - center[1]

  let p = new Array()
  p.push(point)
  let rpoint = multiplyMatrices(p, rotationMatrix(theta))[0]
  // Transpose back
  rpoint[0] = rpoint[0] + center[0]
  rpoint[1] = rpoint[1] + center[1]
  return rpoint
}

export let Point = function(x, y, z) {
    this.x = x
    this.y = y
    this.z = z

    //  Get functions
    this.clone = function() {
      return new Point(this.x, this.y, this.z)
    }
    this.clone_noisy = function(n) {
      return new Point(this.x+Rz(n), this.y+Rz(n), this.z+Rz(n))
    }
    this.getVect = function(i) {
      if (i==2) {
        return [this.x, this.y]
      } else {
        return [this.x, this.y, this.z]
      }
    }

    // Algebra
    this.mult = function(mulVal) {
      this.x = this.x * mulVal
      this.y = this.y * mulVal
    }
    this.div = function(divVal) {
      this.x = this.x / divVal
      this.y = this.y / divVal
    }
    this.div3 = function(divVal) {
      this.x = this.x / divVal
      this.y = this.y / divVal
      this.z = this.z / divVal
    }

    // Norms
    this.norm = function() {
      this.norm = Math.sqrt(this.x**2 + this.y**2)
      this.x = this.x / this.norm
      this.y = this.y / this.norm
    }
    this.norm3 = function() {
      this.norm = Math.sqrt(this.x**2 + this.y**2 + this.z**2)
      this.x /= this.norm
      this.y /= this.norm
      this.z /= this.norm
    }
}

export function CatmullRomsplines(anchorPoints, tho=1, precision=5, close=false, pointType=true) {
  // Add points as a loop
  // Todo: Find other way to create point in case of straight fwd line
  if (close) {
    let lastOne = anchorPoints[anchorPoints.length -1]
    anchorPoints.push(anchorPoints[0])
    anchorPoints.push(anchorPoints[1])
    anchorPoints.unshift(lastOne)
  }

  // Get the CatmullRom Matrix according to 'tho'
  let baseMat = CatmullRom(tho)

  // Filling the points in between ancord points
  let tracePoints = []
  for (var i = 1; i < anchorPoints.length-2; i++) {
    // Points around the section i and i+1
    let points = [
      anchorPoints[i-1],
      anchorPoints[i],
      anchorPoints[i+1],
      anchorPoints[i+2],
    ]
    if (pointType) {
      points = [
        anchorPoints[i-1].getVect(2),
        anchorPoints[i].getVect(2),
        anchorPoints[i+1].getVect(2),
        anchorPoints[i+2].getVect(2),
      ]
    }

    // Todo: Use the distance between two points to create the precision
    // let dist = distance(points[1], points[2])
    for (var j = 0; j<precision; j++){
      let t = j/precision
      let tDer = [1, t, t**2, t**3]
      let calcul = multiplyMatrices(baseMat, points)

      calcul = multiplyMatrices([tDer], calcul)
      if (pointType) {
        let np = new Point(calcul[0][0], calcul[0][1], 0)
        tracePoints.push(np)
      } else {
        tracePoints.push([calcul[0][0], calcul[0][1]])
      }
    }
  }

  // Close the loop, or maye don't
  if (close) {
    tracePoints.push(tracePoints[0])
    tracePoints.push(tracePoints[1])
    tracePoints.unshift(tracePoints[tracePoints.length - 3])
  } else {
    // tracePoints.push(tracePoints[0])
    // tracePoints.push(tracePoints[1])
    // tracePoints.unshift(tracePoints[tracePoints.length - 2])
  }
  return tracePoints;
}

export let setCapital = function(str) {
  return str.charAt(0).toUpperCase() + str.slice(1)
}

export function createTriangulation(p, transp) {
      // Initialize an array to store the points for the Voronoi diagram
    let points = [];

    // Set the number of points to generate
    let numPoints = 100;

    // Generate the points for the Voronoi diagram
    for (let i = 0; i < numPoints; i++) {
      points.push([R(p.width), R(p.height)]);
    }

    // Create the Voronoi diagram
    let voronoi = d3.voronoi()
      .extent([[0, 0], [p.width, p.height]])
      .x(d => d[0])
      .y(d => d[1]);
    let diagram = voronoi(points);
    diagram.polygons().forEach(polygon => {
      var contour = []
      for (let i=0; i< polygon.length; i++) {
        contour.push(new poly2tri.Point(polygon[i][0], polygon[i][1]))
      }
      var swctx = new poly2tri.SweepContext(contour);
      p.push()
      p.beginShape();
      p.fill(150, R()*transp)
      // p.noStroke()
      p.stroke(0)
      p.strokeWeight(0)
      // p.fill(150)

      swctx.triangulate();
      let triangles = swctx.getTriangles();
      // console.log(triangles)
      triangles.forEach(triangle => {
        triangle.getPoints().forEach(point => {
          p.vertex(point.x, point.y);
        });
      });
      p.endShape();
      p.pop()
    });
}

export function drawPencil(points, nAmt, intensity, p) {
    for (let i=1; i<points.length; i++) {
        drawLinePencil([points[i-1], points[i]], nAmt, i, points.length,intensity, p)
    }
}

export function drawLinePencil(points, nAmt, i, l, intensity, p) {
  p.push()
  // p.strokeWeight(1)
  // p.noStroke()
  // p.stroke(87, 77, 68, 50)
  // p.fill(87, 77, 68, 50)
  // p.line(points[0].x, points[0].y, points[1].x, points[1].y)
  let dist = distanceP(points[0], points[1])
  let wid = nAmt*i/l
  let nbP = intensity*5*dist

  for (let t=0; t<nbP; t++) {
    let center = Mp(t/nbP, points[0], points[1])
    let theta = R()*Math.PI*2
    let newP = new Point(
      center.x + Math.cos(theta)*wid,
      center.y + Math.sin(theta)*wid)
    p.circle(newP.x, newP.y, 0.5)
    // p.line(newP.x, newP.y,center.x, center.y, 0.1)
  }
  p.pop()
}


export let distance = function([x, y], [xx, yy]){
  // Euclidian distance in 2d world
  return Math.sqrt(Math.pow(x - xx, 2) +  Math.pow(y - yy, 2));
}
export let distanceP = function(p1, p2){
  // Euclidian distance in 2d world
  return Math.sqrt(Math.pow(p1.x - p2.x, 2) +  Math.pow(p1.y - p2.y, 2));
}
export let subP = function(p1, p2){
  // Euclidian distance in 2d world
  return new Point(p1.x - p2.x, p1.y - p2.y)
   // Math.sqrt(Math.pow(p1.x - p2.x, 2) +  Math.pow(p1.y - p2.y, 2));
}
export let divP = function(p1, value) {
  return new Point(p1.x/value, p1.y/value)
}
export let multP = function(p1, value) {
  return new Point(p1.x*value, p1.y*value)
}
export let addP = function(p1, p2) {
  return new Point(p1.x + p2.x, p1.y + p2.y)
}
export function norm2dP(vector) {
  return Math.sqrt(vector.x**2 + vector.y**2)
}

// Perlin noise map
export function Perlin(nbNodes, size, nbLevels) {
  this.nbLevels = nbLevels
  this.nbNodes = nbNodes
  this.size = size
  this.rand_vect = function(){
      let theta = Math.random() * 2 * Math.PI;
      return {x: Math.cos(theta), y: Math.sin(theta)};
  }
  this.dot_prod_grid = function(x, y, vx, vy){
      let g_vect;
      let d_vect = {x: x - vx, y: y - vy};
      if (this.gradients[[vx,vy]]){
          g_vect = this.gradients[[vx,vy]];
      } else {
          g_vect = this.rand_vect();
          this.gradients[[vx, vy]] = g_vect;
      }
      return d_vect.x * g_vect.x + d_vect.y * g_vect.y;
  }
  this.smootherstep = function(x){
      return 6*x**5 - 15*x**4 + 10*x**3;
  }
  this.interp = function(x, a, b){
      return a + this.smootherstep(x) * (b-a);
  }

  this.gradients = {};
  this.memory = {};

  this.get = function(x, y, ti) {
    // create strates of noise
    let v = 0
    for (let l=1; l<this.nbLevels+1; l++) {
      // if (l==this.nbLevels) {
      //   v = v + this.getLevel(x, y, ti, l)*(1/l)
      // } else {
      v = v + this.getLevel(x, y, ti, l)*(1/(l+1))
      // }
    }
    return v
  }

  this.getLevel = function(x, y, ti, level) {

    x = this.nbNodes*level* (x+level*this.size)/(this.size) + ti/10000
    y = this.nbNodes*level* (y+level*this.size)/(this.size) + ti/10000
    if (this.memory.hasOwnProperty([x,y]))
        return this.memory[[x,y]];
    let xf = Math.floor(x);
    let yf = Math.floor(y);
    //interpolate
    let tl = this.dot_prod_grid(x, y, xf,   yf);
    let tr = this.dot_prod_grid(x, y, xf+1, yf);
    let bl = this.dot_prod_grid(x, y, xf,   yf+1);
    let br = this.dot_prod_grid(x, y, xf+1, yf+1);
    let xt = this.interp(x-xf, tl, tr);
    let xb = this.interp(x-xf, bl, br);
    let v = this.interp(y-yf, xt, xb);
    this.memory[[x,y]] = v;
    return v;
  }
}
