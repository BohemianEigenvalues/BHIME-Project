# BHIME-Project

BHIME (pronounced Bohemian) is an acronym for Bounded Height Integer Matrix Eigenvalues. This project originated in exploring and visualizing the distributions of bounded height integer matrices. Since then, the project has evolved to include matrices that are not necessarily integer matrices. The name has stuck and we typically call the eigenvalues in these images Bohemian Eigenvalues.

There are 4 main parts to the code for this project:
- Generating random matrices
- Sampling eigenvalues from classes of random matrices
- Processing the eigenvalues into a grid in the complex plane
- Producing an image of the eigenvalues

## Generating Random Matrices
Several functions have been provided that will return random matrices given some input values (size, sampling values, etc.). All functions for generating random matrices can be found in the `matrixGenerators` directory. A few simple functions for generating random matrices have been provided.

### `randomMatrix`
The `randomMatrix` function in the `matrixGenerators` directory generates random matrices of a given size with elements sampled uniformly from a given vector of possible values.

__Example__
```matlab
n = 4;  % 4x4 matrix

% Entries for matrix
population = [-1, 0, 1];

% Returns a random 4x4 matrix with entries sampled from [-1, 0, 1];
randomMatrix(population, n);
```

### `randomSymmetricMatrix`
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