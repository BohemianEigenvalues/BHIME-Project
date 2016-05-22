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

__[Full examples](https://github.com/steventhornton/BHIME-Project#examples) as well as the images they produce can be found in the last section__

#### Requirements
- The __Parallel Computing Toolbox__ must be installed. If you do not have the parallel computing toolbox install you may modify the code slightly by changing any `parfor` loops to `for` loops.
- All contents (files and subfolders) of the `src` directory must be in your Matlab working directory.

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
- It must return square matrices of the same size each time you call the generator


# Sampling Eigenvalues

There is one main function for generating sample eigenvalues, the `generateRandomSample` function in the `dataGeneration` directory.
This function will generate a large random sample of matrices using the provided generating function, compute their eigenvalues and store them in `.mat` files.

This function will:
- Create a new directory `Data` in the input working directory
- By default it will generate 1000000/n (n is the size of your matrices) matrices and store their eigenvalues and condition numbers in a `.mat` file

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

The function used to process the eigenvalues into a grid is the `processData` function in the `dataProcessing` directory.

This function will:
- Create a new directory `ProcessedData` in the input working directory

## Examples

#### Simple example
```matlab
workingDir = '~/ComplexSymmetric/';

% Process by density
fname = processData(workingDir);

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

The function for converting the processed data into an image is the `processImage` function in the `imageProcessing` directory.

This function will:
- Create a new directory `Images` in the input working directory
- Write an image (.png) in the Images directory
- Automatically name the image to avoid overwriting any other images

# Options
Each of the methods

- `generateRandomSample`
- `processData`
- `processImage`

take an optional input value. It is a Matlab struct that controls several things when producing/processing eigenvalues. All options are summarized in the table below.

| Option Name | Default | Details |
| ----------- | ------- | ------- |
| `backgroundColor` | `[0, 0, 0]` (black) | Set this to a vector with 3 values representing the RGB values (between 0 and 1) to change the background color of the image that is produced |
| `colorByCond` | `false` | When set to `true`, the colors in the output image represent the average condition number in that pixel. To use this option the condition number data must be computed. Use the `computeCond` option with the `generateRandomSample` or `generateAllMatrices` functions to compute the eigenvalues condition numbers. |
| `computeCond` | `false` | When set to `true`, the `generateRandomSample` or `generateAllMatrices` function will compute and store the condition numbers of the eigenvalues. The functions will take roughly 6 times longer to execute when set to `true`. |
| `dataPrecision` | `'single'` | This can be set to either `'single'` or `'double'`. This option __does not__ affect the precision the eigenvalues are computed in. Eigenvalues are __always__ computed in double precision. By setting this to `'single'` the eigenvalues will be __stored__ in single precision. That is, they will be computed in double precision and cast to single precision for storage. |
| `filenamePrefix` | `'BHIME'` | The name that will be used when naming the data files. The names of the data files take the form: `filenamePrefix + '_' + i` where `i` is a positive integer. |
| `height` | 1001 (pixels) | The height (in pixels) of the image to be produced. The width is determined from the `margin` such that each grid point is square. |
| `ignoreReal` | `false` | Set this value to `true` to ignore any eigenvalues where the imaginary part of the eigenvalue is within the `ignoreRealTol` of zero. This option is used for the `processData` function. |
| `ignoreRealData` | `false` | When this option is `true`, eigenvalues where the imaginary part is within `ignoreRealTol` of zero will not be stored in the data file. This option is used for `generateRandomSample` and `generateAllMatrices`. By setting this to `true` you can substantially reduce the size of the data files. |
| `ignoreRealTol` | `1e-10` | A complex number `z` will be considered a real value, and therefore ignored (if `ignoreReal` or `ignoreRealData` is `true`) if `abs(Im(z)) < ignoreRealTol`. |
| `map` | `@(z) z` (no mapping) | Map the eigenvalues by a given function handle. __Must be vectorized.__ |
| `margin` | Large enough to fit all the data points in the first data file. | Must be a struct with keys: <ul><li>`bottom`</li><li>`top`</li><li>`left`</li><li>`right`</li></ul> that indicate the margins for the image. |
| `matricesPerFile` | `floor(1e6/matrixSize)` | Control how many matrices eigenvalues (and condition numbers) are in each data file. |
| `maxDensity` | 0 (i.e. not set) | Controls the maximum eigenvalue count used for coloring the image. Any point with a higher density than this this value is colored as the final color in the colormap. |
| `minDensity` | 0 (i.e. not set) | Controls the minimum eigenvalue count used for coloring the image. Any point with a lower density than this this value is colored as the final color in the colormap. |
| `numCharPolyFiles` | 1 | The number of characteristic polynomial (plain text files) to convert to eigenvalue data files. |
| `numDataFiles` | 1 | Set this option to a positive integer if you would like to generate multiple files with data where each file contains the eigenvalues and their condition numbers for `matricesPerFile` random matrices. |
| `numProcessFiles` | Number of files in the `Data` directory | The number of data files for the `processData` function to use. |
| `outputFileType` | `'mat'` | Can set to `'txt'` if you want the processed data written to a text file. |
| `overrideDataDir` | `''` | Set this option to a non-empty string indicating a directory to write the data files to. If not set they will be written to the directory `workingDir/Data/`. If the directory does not exist the `generateRandomSample` and `generateAllMatrices` functions will create the directory. |
| `overrideImagesDir` | `''` | Set this option to a non-empty string indicating a directory to write the images to. If not set they will be written to the directory `workingDir/Images/`. If the directory does not exist the `processImage` function will create it. |
| `overrideProcessedDataDir` | `''` | Set this option to a non-empty string indicating a directory to write the processed data files to. If not set they will be written to the directory `workingDir/ProcessedData/`. If the directory does not exist the `processData` function will create it. |
| `startFileIndex` | 1 greater than the  highest index of the files in the data directory or 1 if no files have been written | Only use this if you have already computed data and would like to compute more. |
| `storeDataWithSymmetry` | `false` | When this option is set to `true`, any values not in the upper right quadrant of the complex plane will not be stored in the output data files. That is, only values where `Im(z) >= 0` and `Re(z) >= 0` are stored. This option is used by `generateRandomSample` and `generateAllMatrices`. |
| `symmetryIm` | `false` | If `true`, symmetry across the imaginary axis will be used to effectively double the number of points. |
| `symmetryRe` | `false` | If `true`, symmetry across the real axis will be used to effectively double the number of points. |

__How to determine a good value for `matricesPerFile`__:
Each file will use `32*matrixSize*matricesPerFile` bits if the condition numbers are not computed, otherwise each file will be `64*matrixSize*matricesPerFile` bits, make sure this value is less than the amount of RAM your computer has. Ideally this value is less than `1/numCores` where `numCores` is the number of processors you are allowing the parallel computing toolbox in Matlab to use.
By setting the `numDataFiles` option you can generate many files, each of which will contain data on `matricesPerFile` random matrices.
With 16GB of RAM, the maximum value I use is `matricesPerFile = 1e7`.

# Examples

### Example 1
This simple example explores the eigenvalues of random 5x5 matrices with entries sampled from {-1, 0, 1}
```matlab
workingDir = '~/Real5x5/';

% Generate Data ------------------

% The generator (5x5 matrices with entries sampled from {-1, 0, 1})
g = @() randomMatrix([-1, 0, 1], 5);

% Generate the random sample
generateRandomSample(g, workingDir);

% Process Data -------------------

pFilename = processData(workingDir);

% Make the image ------------------

% Colormap
T = [  0,   0,   0;
       0,   0, 255;
       0, 255, 255;
       0, 255,   0;
     255, 255,   0;
     255,   0,   0;
     255,   0,   0;
     255, 255, 255;
     255, 255, 255]/255;

% Weights for the colormap
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% Plot image ---------------------
plotImage(workingDir, pFilename,  T, x);
```

This produces the image:

<p align="center">
    <img alt="5x5 matrices with entries sampled from {-1, 0, 1}" src="https://s3.amazonaws.com/stevenethornton.github/Real5x5.png"/>
</p>

### Example 2
```matlab
n = 4;  % 4x4 matrices

% Entries for matrices
population = exp(2i*pi/5*(0:4));

% The generator
g = @() randomSymmetricMatrix(population, n);

workingDir = '~/ComplexSymmetric/';

% Set the options
margin = struct('bottom', -4, ...
                   'top',  4, ...
                  'left', -4, ...
                 'right',  4);
opts = struct('numDataFiles', 10, ...
           'matricesPerFile', 1e6, ...
                    'height', 501, ...
                    'margin', margin);

% Generate the data (may take a few minutes)
generateRandomSample(g, workingDir, opts);

% Process the data
fname = processData(workingDir, opts);

% Colormap
T = [  0,   0,   0;
       0,   0, 255;
       0, 255, 255;
       0, 255,   0;
     255, 255,   0;
     255,   0,   0;
     255,   0,   0;
     255, 255, 255;
     255, 255, 255]/255;

% Weights for colormap
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% Make an image
processImage(workingDir, fname, T, x);
```
produces the image:

<p align="center">
    <img alt="4x4 Complex Symmetric Matrices" src="https://s3.amazonaws.com/stevenethornton.github/ComplexSymmetric_4x4.png"/>
</p>

### Example 3
```matlab
n = 5;  % 5x5 matrices

% Entries for matrices
population = [-20, -1, 0, 1, 20];

% The generator
g = @() randomMatrix(population, n);

workingDir = '~/Real5x5-wide/';

% Set the options
margin = struct('bottom', -40, ...
                   'top',  40, ...
                  'left', -50, ...
                 'right',  50);
opts = struct('matricesPerFile', 1e6, ...
                       'height', 501, ...
                       'margin', margin);

% Generate the data (may take a few minutes)
generateRandomSample(g, workingDir, opts);

% Process the data
fname = processData(workingDir, opts);

% Colormap
T = [  0,   0,   0;
       0,   0, 255;
       0, 255, 255;
       0, 255,   0;
     255, 255,   0;
     255,   0,   0;
     255,   0,   0;
     255, 255, 255;
     255, 255, 255]/255;

% Weights for colormap
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% Make an image
processImage(workingDir, fname, T, x);
```
produces the image:

<p align="center">
    <img alt="5x5 Symmetric Complex Matrices" src="https://s3.amazonaws.com/stevenethornton.github/Real5x5-Wide.png"/>
</p>

### Example 4
```matlab
A = [0,  0,  0,  0;
    -1, -1,  1,  0;
     0,  0,  0,  0;
    -1, -1, -1, -1];

% The entries that will follow a uniform distribution A(1,4) and A(3,1)
entries = [1, 4;
           3, 1];

% The generator
g = @() givenMatrixUniform(A, entries, -5, 5);

workingDir = '~/Given_4x4/';

% Set the options
margin = struct('bottom', -4, ...
                   'top',  4, ...
                  'left', -4, ...
                 'right',  4);
opts = struct('matricesPerFile', 1e6, ...
                       'height', 501, ...
                       'margin', margin, ...
              'backgroundColor', [1, 1, 1]);

% Generate the data (may take a few minutes)
generateRandomSample(g, workingDir, opts);

% Process the data
fname = processData(workingDir, opts);

% Purple gradient
T = [225, 128, 200;
       0,   0,   0]/255;
x = [0.0, 1.0];

% Make an image
processImage(workingDir, fname, T, x, opts);
```
produces the image:

<p align="center">
    <img alt="4x4 matrix with 2 uniformly sampled entries." src="https://s3.amazonaws.com/stevenethornton.github/Given_4x4.png"/>
</p>



# To Do

- Add method for automatic colormap weights
- Get the condition number coloring in `processData` working