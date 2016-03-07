# BHIME-Project

BHIME (pronounced Bohemian) is an acronym for Bounded Height Integer Matrix Eigenvalues. This project originated in exploring and visualizing the distributions of bounded height integer matrices. Since then, the project has evolved to include matrices that are not necessarily integer matrices. The name has stuck and we typically call the eigenvalues in these images Bohemian Eigenvalues.

There are 4 main parts to the code for this project:
- Generating random matrices
- Sampling eigenvalues from classes of random matrices
- Processing the eigenvalues into a grid in the complex plane
- Producing an image of the eigenvalues

## Generating Random Matrices
Several functions have been provided that will return random matrices given some input values (size, sampling values, etc.). All functions for generating random matrices can be found in the `matrixGenerators` directory. A few simple functions for generating random matrices have been provided.

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
This function will generate a large random sample of matrices using the provided generating function, compute their eigenvalues and store them in files.

This function will:
- Create a new directory `Data` in the input working directory
- By default it will generate 1000000/n (n is the size of your matrices) matrices and store their eigenvalues and condition numbers in a `.mat` file
- 

### Options
There are 4 options that can be provided to the function.

#### `filenamePrefix`
This option is set to `'BHIME'` by default.

It is the name that will be used when naming the data files. The names of the data files take the form: `filenamePrefix + '_' + i` where `i` is a positive integer.

#### `startFileIndex`
This option is set to 1 more than the highest index of the files in the data directory

#### `numFiles`
This is set to 1 by default

Set this option to a positive integer if you would like to generate multiple files with data




Example of generating eigenvalue data
```matlab
n = 4;  % 4x4 matrices

% Entries for matrices
population = exp(2i*pi/5*(0:4));

% The generator
g = @() randomSymmetricMatrix(population, n);

workingDir = '~/ComplexSymmetric/';

generateRandomSample(g, workingDir);
```