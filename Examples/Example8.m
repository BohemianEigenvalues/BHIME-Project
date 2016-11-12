% ----------------------------------------------------------------------- %
% EXAMPLE 8 -- Eigenvalues of 20x20 companion matrices of polynomials     %
%              with the coefficients sampled uniformly from the set       %
%              {-1, 1}.                                                   %
% ----------------------------------------------------------------------- %

% Matrix size (polyDeg = 20)
n = 19;

% Entries for matrices
population = [-1, 1];

% The generator
g = @() randomCompanionMatrix(population, n);

workingDir = '~/CompanionMatrix/';

% Set the options
margin = struct('bottom', -1.6, ...
                   'top',  1.6, ...
                  'left', -1.6, ...
                 'right',  1.6);

opts = struct('matricesPerFile', 1e6, ...
                 'numDataFiles', 1, ...
                       'height', 2001, ...
                       'margin', margin);

% Generate the data (may take a few minutes)
generateAllMatrices(g, workingDir, length(population)^n, opts);

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

% The weights for the colormap
x = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0];

% Make an image
processImage(workingDir, fname, T, x);