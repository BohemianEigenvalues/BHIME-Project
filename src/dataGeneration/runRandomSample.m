% runRandomSample

% Matrix Size
n = 4;

% Entries for matrix
% population = [-5i, -1, 0, 1, 5i];

% population = exp(1i*pi*(1:2:7)/4);

% population = exp(1i*pi*(2/3)*[0,1,2]);
population = exp(1i*pi/7 * (0:2:13));

% The generator
g = @() randomSymmetricMatrix(population, n);

workingDir = '~/Desktop/TestCode7/';

options = struct('numFiles',10);

generateRandomSample(g, workingDir, options);