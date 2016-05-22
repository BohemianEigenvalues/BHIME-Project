# BHIME Project

BHIME (pronounced Bohemian) is an acronym for Bounded Height Integer Matrix Eigenvalues. 
This project originated in exploring and visualizing the distributions of the eigenvalues of bounded height integer matrices. 
The project has since evolved to include many other classes of matrices.
The name has stuck and we typically call the eigenvalues in these images *Bohemian Eigenvalues*.

Slides on the Bohemian Eigenvalue Project can be found [here](https://s3.amazonaws.com/stevenethornton.github/BHIME+Slides.pdf).

__Full examples as well as the images they produce can be found in the [Examples section](https://github.com/steventhornton/BHIME-Project#examples) or in the [Examples folder]().__

#### Matlab Requirements
- The __Parallel Computing Toolbox__ must be installed. If you do not have the parallel computing toolbox install you may modify the code slightly by changing all `parfor` loops to `for` loops.
- All contents (files and subfolders) of the `src` directory must be in your Matlab working directory.

Detailed instructions on how the code works can be found in the [readme]() in the `src` directory, here I will only give examples of its use as well as the full list of options.

# Examples

## Example 1
The eigenvalues of a random sample of 5x5 matrices where the entries are sampled uniformly from {-1, 0, 1}.
See [`Examples/Example1.m`]() for a more detailed explanation of how this example works.
```matlab
% Set the working directory
workingDir = '~/Real5x5_d3/';

% The generator (5x5 matrices with entries sampled from {-1, 0, 1})
g = @() randomMatrix([-1, 0, 1], 5);

% Generate a random sample
generateRandomSample(g, workingDir);

% Process the sample into a grid in the complex plane
pFilename = processData(workingDir);

% The colormap for the resulting image
T = [  0,   0,   0;
       0,   0, 255;
       0, 255, 255;
       0, 255,   0;
     255, 255,   0;
     255,   0,   0;
     255,   0,   0;
     255, 255, 255;
     255, 255, 255]/255;

% The weights for the colormap
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% Create the image
processImage(workingDir, pFilename, T, x);
```

This produces the image:

<p align="center">
    <img alt="5x5 matrices with entries sampled from {-1, 0, 1}" src="https://s3.amazonaws.com/stevenethornton.github/Real5x5.png"/>
</p>

# Options
Each of the methods

- `generateRandomSample`
- `generateAllMatrices`
- `processCharPolyFile`
- `processData`
- `processImage`

take an optional final input value. It is a Matlab struct that controls several things when producing/processing eigenvalues. All options are summarized in the table below.

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