# BHIME Project

BHIME (pronounced Bohemian) is an acronym for Bounded Height Integer Matrix Eigenvalues. 
This project originated in exploring and visualizing the distributions of the eigenvalues of bounded height integer matrices. 
The project has since evolved to include many other classes of matrices.
The name has stuck and we typically call the eigenvalues in these images *Bohemian Eigenvalues*.

The code is separated into 4 main components:
- [Generating random matrices](https://github.com/steventhornton/BHIME-Project#generating-random-matrices)
- [Sampling eigenvalues from classes of random matrices](https://github.com/steventhornton/BHIME-Project#sampling-eigenvalues)
- [Processing the eigenvalues into a grid in the complex plane](https://github.com/steventhornton/BHIME-Project#processing-the-eigenvalues)
- [Producing an image of the eigenvalues](https://github.com/steventhornton/BHIME-Project#making-an-image)

# Generating Random Matrices
Several functions have been provided that will return random matrices given some input values (size, sampling values, etc.). All functions for generating random matrices can be found in the `matrixGenerators` directory.

### Random Matrices
The `randomMatrix` function in the `matrixGenerators` directory generates random matrices of a given size with elements sampled uniformly from a given vector of possible values.

__Example__
```matlab
n = 4;  % 4x4 matrix

% Entries for matrix
population = [-1, 0, 1];

% Returns a random 4x4 matrix with entries sampled from [-1, 0, 1];
randomMatrix(population, n);
```

### Random Symmetric Matrices
The `randomSymmetricMatrix` function in the `matrixGenerators` directory generates random symmetric matrices of a given size with elements sampled uniformly from a given vector of possible values.

__Example__
```matlab
n = 4;  % 4x4 matrix

% Entries for matrix
population = [-1i, 0, sqrt(2)];

% Returns a random symmetric 4x4 matrix with entries sampled from [-i, 0, sqrt(2)];
randomSymmetricMatrix(population, n);
```

### Custom Generator Functions
Writing your own random matrix generator function is simple. There is only one rule:
- They must return square matrices of the same size each time you call the generator


# Sampling Eigenvalues

There is one main function for generating sample eigenvalues, the `generateRandomSample` function in the `dataGeneration` directory.
This function will generate a large random sample of matrices using the provided generating function, compute their eigenvalues and store them in `.mat` files.

This function will:
- Create a new directory `Data` in the input working directory
- By default it will generate 1000000/n (n is the size of your matrices) matrices and store their eigenvalues and condition numbers in a `.mat` file

### Options
There are 4 options that can be provided to the function.

| Option Name | Default | Details |
| ----------- | ------- | ------- |
| `filenamePrefix` | `'BHIME'` | The name that will be used when naming the data files. THe names of the data files take the form: `filenamPrefix + '_' + i` where `i` is a positive integer. |
| `startFileIndex` | 1 more than the highest index of the files in the data  directory, 1 if no files have been written | Only use this if you have already computed data and would like to compute more |
| `numFiles` | 1 | Set this option to a positive integer if you would like to generate multiple files with data where each file contains the eigenvalues and their condition numbers for `matricesPerFile` random matrices |
| `matricesPerFile` | `1000000/matrixSize` | Control how many matrices eigenvalues/condition numbers are in each file |

__How to determine a good value for `matricesPerFile`__:
Each file will use `64*matrixSize*matricesPerFile` bits, make sure this value is less than the amount of RAM your computer has.
By setting the `numFiles` option you can generate many files, each of which will contain data on `matricesPerFile` random matrices.

### Examples

__Simple example of generating eigenvalue data__
```matlab
n = 4;  % 4x4 matrices

% Entries for matrices
population = exp(2i*pi/5*(0:4));

% The generator
g = @() randomSymmetricMatrix(population, n);

workingDir = '~/ComplexSymmetric/';

generateRandomSample(g, workingDir);
```

__Example with options__
```matlab
n = 5;  % 5x5 matrices

% Entries for matrices
population = [-1, 0, 1];

% The generator
g = @() randomSymmetricMatrix(population, n);

workingDir = '~/Real5x5/';

% Options
options = struct('filenamePrefix', 'foobar', ...
                 'startFileIndex', 1, ...
                 'numFiles', 10, ...
                 'matricesPerFile', 1e5);

generateRandomSample(g, workingDir, options);
```

# Processing the Eigenvalues

## Options

| Option Name | Default | Details |
| ----------- | ------- | ------- |
| `height` | 1001 (pixels) | The height (in pixels) of the grid to be used. The width is determined from the `margin` such that each grid point is square. |
| `margin` | Large enough to fit all the points in the first data file | Must be a struct with keys: <ul><li>`bottom`</li><li>`top`</li><li>`left`</li><li>`right`</li></ul> that indicate the margins for the image. |
| `dataFilePrefix` | `BHIME` | The prefix for the .mat files that contain the eigenvalues and their condition numbers |
| `outputFileType` | `mat` | Can set to `txt` if you want the processed data written to a text file |
| `symmetry` | `false` | If `true`, symmetry across the real and imaginary axes will be used to effectively quadruple the number of points |
| `numFiles` | All files in the `Data` directory | The number of data files to process |
| `map` | `@(z) z` (no mapping) | Map the eigenvalues by a given function handle. __Must be vectorized.__ |

## Examples

#### Simple example
```matlab
workingDir = '~/ComplexSymmetric/';

% Process by density
fname = processData(workingDir, 'density');

% fname will be the name of the file containing the data
```

#### Example with options
```matlab
workingDir = '~/Real5x5/';

margin = struct('bottom', -4, ...
                'top',     4, ...
                'left',   -4, ...
                'right',   4);

opts = struct('height', 2001, ...
              'margin', margin, ...
              'dataFilePrefix', 'foobar', ...
              'outputFileType', 'txt', ...
              'symmetry', true, ...
              'numFiles', 5, ...
              'map', @(z) 1./z);

% Process by condition number 
fname = processData(workingDir, 'cond', opts);
```

# Making an Image

# Complete Examples

### Example 1:
This simple example explores the eigenvalues of random 5x5 matrices with entries sampled from {-1, 0, 1}
```matlab
workingDir = '~/Real5x5/';

% Generate Data ------------------

% The generator (5x5 matrices with entries sampled from {-1, 0, 1})
g = @() randomMatrix([-1, 0, 1], 5);

% Generate the random sample
generateRandomSample(g, workingDir);

% Process Data -------------------

% Color the image by the density of eigenvalues
colorBy = 'density'

pFilename = processData(workingDir, colorBy);

% Make the image ------------------

T = [0, 0,   0,    0,    255,  255,  255, 255,  255;
     0, 0,   255,  255,  255,  0,    0,   255,  255;
     0, 255, 255,  0,    0,    0,    0,   255,  255];
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% Plot image ---------------------
plotImage(workingDir, pFilename,  T, x);
```

This produces the image:
![5x5 Integer Matrices](https://s3.amazonaws.com/stevenethornton.github/Image-1.png)

### Example 2
```matlab
n = 4;  % 4x4 matrices

% Entries for matrices
population = exp(2i*pi/5*(0:4));

% The generator
g = @() randomSymmetricMatrix(population, n);

workingDir = '~/ComplexSymmetric/';

% Set the options
opts = struct('numFiles', 10, ...
              'matricesPerFile', 1e6);

% Generate the data (may take a few minutes)
generateRandomSample(g, workingDir, opts);

% Options for processing the data
margin = struct('bottom', -4, ...
                'top', 4, ...
                'left', -4, ...
                'right', 4);

opts = struct('height', 501, ...
              'margin', margin);

% Process the data
colorBy = 'density';
fname = processData(workingDir, colorBy, opts);

T = [0, 0,   0,    0,    255,  255,  255, 255,  255;
      0, 0,   255,  255,  255,  0,    0,   255,  255;
      0, 255, 255,  0,    0,    0,    0,   255,  255]'./255;
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% Make an image
processImage(workingDir, fname, T, x);
```
produces the image:

<p align="center">
![4x4 Symmetric Complex Matrices](https://s3.amazonaws.com/stevenethornton.github/ComplexSymmetric_4x4.png)
</p>
