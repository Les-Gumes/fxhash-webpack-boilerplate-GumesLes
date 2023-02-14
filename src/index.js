import p5 from "p5";

// these are the variables you can use as inputs to your algorithms
console.log(fxhash)   // the 64 chars hex number fed to your algorithm
console.log(fxrand()) // deterministic PRNG function, use it instead of Math.random()
const seed = ~~(fxrand()*123456789);
let s;

const numCircles = ~~(fxrand()*500) + 100;

window.$fxhashFeatures = {
  "Density": numCircles > 500?"High":(numCircles<200?"Low":"Medium")
}


let sketch = function(p5) {

  p5.setup = function() {
    p5.noLoop();
    s = p5.min(p5.windowWidth, p5.windowHeight);
    p5.createCanvas(s, s);
  };

  p5.draw = function() {
    p5.randomSeed(seed);
    p5.background("#FFD");

  };

  p5.windowResized = function() {
    s = p5.min(p5.windowWidth, p5.windowHeight);
    p5.resizeCanvas(s, s);
  }
}

let myp5 = new p5(sketch, window.document.body);
