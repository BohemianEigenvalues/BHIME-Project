# BHIME-Project

BHIME (pronounced Bohemian) is an acronym for Bounded Height Integer Matrix Eigenvalues. This project originated in exploring and visualizing the distributions of bounded height integer matrices. Since then the project has evolved to include matrices that are not necessarily integer matrices. The name has stuck and we typically call the eigenvalues in these images Bohemian Eigenvalues.

There are 4 main parts to the code for this project
- Generating random matrices
- Sampling eigenvalues from classes of random matrices
- Processing the eigenvalues into a grid in the complex plane
- Producing an image of the eigenvalues


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