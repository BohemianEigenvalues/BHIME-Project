# BHIME-Project

BHIME (pronounced Bohemian) is an acronym for Bounded Height Integer Matrix Eigenvalues. This project originated in exploring and visualizing the distributions of bounded height integer matrices. Since then, the project has evolved to include matrices that are not necessarily integer matrices. The name has stuck and we typically call the eigenvalues in these images Bohemian Eigenvalues.

There are 4 main parts to the code for this project:
- Generating random matrices
- Sampling eigenvalues from classes of random matrices
- Processing the eigenvalues into a grid in the complex plane
- Producing an image of the eigenvalues

## Generating Random Matrices
Several functions have been provided that will return random matrices given some input values (size, sampling values, etc.). The `matrixGenerators` directory contains several random matrix generators that can be used, custom generators can be written (see below).

### Random Matrices
The function `randomMatrix` in the `matrixGenerators` directory can be used to generate random matrices where the elements are sampled uniformly from a given vector of possible values. The function takes two input values, a vector of values to sample from and a positive integer for the size of the matrix.

__Example__
```matlab
n = 4;  % 4x4 matrix

% Entries for matrix
population = [-1, 0, 1];

% Returns a random 4x4 matrix with entries sampled from [-1, 0, 1];
randomMatrix(population, n);
```

### Random Symmetric Matrices
The function `randomSymmetricMatrix` in the `matrixGenerators` directory generates random symmetric matrices given the size of the matrix and a vector to sample values from.

__Example__
```matlab
n = 4;  % 4x4 matrix

% Entries for matrix
population = [-1, 0, 1];

% Returns a random symmetric 4x4 matrix with entries sampled from [-1, 0, 1];
randomSymmetricMatrix(population, n);
```

### Custom Generator Functions
Writing your own random matrix generator function is simple. There is only one rule:
- They must return square matrices of the same size each time you call the generator


Example of generating eigenvalue data
```matlab
n = 4;  % 4x4 matrices

% Entries for matrices
population = exp(1i*pi/5*(0:2:8));

% The generator
g = @() randomSymmetricMatrix(population, n);

workingDir = '~/ComplexSymmetric/';

generateRandomSample(g, workingDir);

```