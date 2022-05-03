// 2D Arrays
var current_grid;
var new_grid;
// Reaction-Diffusion variables
// Diffusion rate of A
var dA = 1.0;
// Diffusion Rate of B
var dB = 0.5;
// Feed rate of A
var f = 0.055;
// Kill rate of B
var k = 0.062;
// Canvas Variables
var canvas_width = 400;
var canvas_height = 400;
var seed_size = 20;

/**
* setup() creates the canvas and the grids that are used to
* hold the concentration values. The grid is initialized such
* that density of Chemical A is 1.0 and Chemical B is 0.0.
* Then, a seed of Chemical B is added.
*/
function setup() {
  // Create a canvas of 200px by 200px
  createCanvas(canvas_width, canvas_height);
  pixelDensity(1);

  // Establish your 2D array variables
  current_grid = [];
  new_grid = [];

  // Loop through the width of the canvas and fill
  // each column with an empty array
  for(var x = 0; x < width; x++){

    // We are doing the same setup for each of our
    // grids to start.
    current_grid[x] = [];
    new_grid[x] = [];

    // For each column, populate the values of the
    // chemical concentrations. The entire grid will be
    // populated with Chemical A to start.
    for(var y = 0; y < height; y++){
      current_grid[x][y] = {
        a: 1,
        b: 0
      }
      new_grid[x][y] = {
        a: 1,
        b: 0
      }
    }
  }

  // Calculate where the B seed should start and end
  // based on the canvas size and the seed size.
  // The B seed is a square.
  var start_x = floor(width/2 - seed_size/2);
  var end_x = start_x + seed_size;
  var start_y = floor(height/2 - seed_size/2);
  var end_y = start_y + seed_size;

  // Loop through the current grid and add the seed.
  for(var x2 = start_x; x2 < end_x; x2++){
    for(var y2 = start_y; y2 < end_y; y2++){
      current_grid[x2][y2].b = 1.0;
    }
  }
}

/**
* This function is the implementation of the Gray-Scott algorithm.
* The current grid is used to populate the new grid.
*/
function updateConcentrations(){
  // Note: We start at 1 and end at length - 1 so that the laplace
  // function does not hit an index out of bounds error
  for (var x = 1; x < width - 1; x++) {
    for (var y = 1; y < height - 1; y++) {
      // If you look at the equations we're using, the first value
      // is the old concentration valaue from the previous timestep.
      var a = current_grid[x][y].a;
      var b = current_grid[x][y].b;

      // Our next values, the diffusion values, we have already set
      // way up top. Remember they are called dA and dB.

      // We'll write our laplace function later, but we'll name them
      // laplaceA(x, y) and laplaceB(x, y).

      // Next, we have our AB^2 term
      var ab2 = a*b*b;

      // This is where things diverge. We'll start first with the
      // feed term for Chemical A
      var feed = f * (1.0-a);

      // Now, let's look at our kill term for Chemical B
      var kill = (k + f)*b;

      // Our timestep will be 1, so no need to use the delta t term.

      // Now, let's update the variables in our "next" grid
      new_grid[x][y].a = a + (dA * laplaceA(x, y) - ab2 + feed);
      new_grid[x][y].b = b + (dB * laplaceB(x, y) + ab2 - kill);
    }
  }
}

/**
* populatePixels() takes the concentration data from the grid
* and transfers it to the pixel grid.
*/
function populatePixels(){
  // https://p5js.org/reference/#/p5/loadPixels
  // for documentation.
  loadPixels();

  // loop through the canvas
  for (var x = 0; x < width; x++) {
    for (var y = 0; y < height; y++) {
      // All of the pixel data is now contained in pixels[].
      // Even though the canvas looks like a 2D Array, the
      // data structure is actually a 1D array. That means we
      // have to take the x and y coordinate and convert it to
      // its corresponding 1D array position using this formula.
      // To add onto this craziness, there are 4 spots for each
      // pixel, containing R, G, B, and A values.
      var pixel_loc = (x + y * width) * 4;

      // Get the new concentration values
      var a = new_grid[x][y].a;
      var b = new_grid[x][y].b;

      // Now we need to actually associate a color with the
      // concentration values. Remember the diagram from
      // above? If A > B, then the pixel will be white. But,
      // if B < A, it will be black. We round the function so
      // the color value is an integer.
      var color = round((a - b) * 255);

      // We now have to make sure we don't have anything above
      // or below 255, as that is not a color value.
       color = constrain(color, 0, 255);

       // We want each pixel to have our color, making the image
       // grayscale. We want opacity to be 100%, so our alpha
       // value will stay at 255.
       pixels[pixel_loc + 0] = color;
       pixels[pixel_loc + 1] = color;
       pixels[pixel_loc + 2] = color;
       pixels[pixel_loc + 3] = 255;
    }
  }

  // https://p5js.org/reference/#/p5/updatePixels
  // for documentation.
  updatePixels();
}

/**
* swapGrids() switches the data contained in current_grid and new_grid
*/
function swapGrids() {
  var temp_grid = current_grid;
  current_grid = new_grid;
  new_grid = temp_grid;
}

/**
* The Laplacian function for Chemical A.
* @param x: The X coordinate for the current cell.
* @param y: The Y coordinate for the current cell.
* @return the combined value of the surrounding cells.
*/
function laplaceA(x, y) {
  var sumA = 0;

  // From before: center weight -1, adjacent neighbors .2, diagonals .05.

  // Center weight -1
  sumA += current_grid[x][y].a * -1;

  //Adjacent neighbors weight = 0.2
  // Left of center
  sumA += current_grid[x - 1][y].a * 0.2;
  // Right of center
  sumA += current_grid[x + 1][y].a * 0.2;
  // Below center
  sumA += current_grid[x][y + 1].a * 0.2;
  // Above center
  sumA += current_grid[x][y - 1].a * 0.2;

  // Diagonal weights .05
  // Above Left
  sumA += current_grid[x - 1][y - 1].a * 0.05;
  // Above Right
  sumA += current_grid[x + 1][y - 1].a * 0.05;
  // Below Right
  sumA += current_grid[x + 1][y + 1].a * 0.05;
  // Below Left
  sumA += current_grid[x - 1][y + 1].a * 0.05;
  return sumA;
}

/**
* The Laplacian function for Chemical B
* @param x: The X coordinate for the current cell.
* @param y: The Y coordinate for the current cell.
* @return the combined value of the surrounding cells.
*/
function laplaceB(x, y) {
  var sumB = 0;

  // From before: center weight -1, adjacent neighbors .2, diagonals .05.

  // Center weight -1
  sumB += current_grid[x][y].b * -1;

  //Adjacent neighbors weight = 0.2
  // Left of center
  sumB += current_grid[x - 1][y].b * 0.2;
  // Right of center
  sumB += current_grid[x + 1][y].b * 0.2;
  // Below center
  sumB += current_grid[x][y + 1].b * 0.2;
  // Above center
  sumB += current_grid[x][y - 1].b * 0.2;

  // Diagonal weights .05
  // Above Left
  sumB += current_grid[x - 1][y - 1].b * 0.05;
  // Above Right
  sumB += current_grid[x + 1][y - 1].b * 0.05;
  // Below Right
  sumB += current_grid[x + 1][y + 1].b * 0.05;
  // Below Left
  sumB += current_grid[x - 1][y + 1].b * 0.05;
  return sumB;
}

/**
* Performs the 3 steps necessary at each timestep.
*/
function draw() {
  updateConcentrations();
  populatePixels();
  swapGrids();
}
