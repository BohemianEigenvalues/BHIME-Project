% ----------------------------------------------------------------------- %
% EXAMPLE 4 -- Eigenvalues of the matrix                                  %
%                 [ 0  0  0  A]                                           %
%                 [-1 -1  1  0]                                           %
%                 [ B  0  0  0]                                           %
%                 [-1 -1 -1 -1]                                           %
%              where A and B are sampled independently from a continuous  %
%              uniform distribution on (-5, 5).                           %
% ----------------------------------------------------------------------- %

workingDir = '~/Given_4x4/';

A = [0,  0,  0,  0;
    -1, -1,  1,  0;
     0,  0,  0,  0;
    -1, -1, -1, -1];

% The entries that will follow a uniform distribution A(1,4) and A(3,1)
entries = [1, 4;
           3, 1];

% The generator
g = @() givenMatrixUniform(A, entries, -5, 5);

% Margin
margin = struct('bottom', -4, ...
                'top',     4, ...
                'left',   -4, ...
                'right',   4);

% The background color of the output image will be blue
opts = struct('matricesPerFile', 1e6, ...
              'height',          501, ...
              'margin',          margin, ...
              'backgroundColor', [28, 38, 137]/255);

% Generate the data (may take a few minutes)
generateRandomSample(g, workingDir, opts);

% Process the data
pFilename = processData(workingDir, opts);

% Gradient
T = [28, 38, 137;
     227, 227, 238;
     255, 255, 255]/255;
x = [0.0, 0.6, 1.0];

% Make an image
processImage(pFilename, workingDir, T, x, opts);