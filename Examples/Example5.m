% ----------------------------------------------------------------------- %
% EXAMPLE 5 -- Eigenvalues of 5x5 matrices with entries in                %
%              {-1, -1/10000, 0, 1/10000, 1}                              %
% ----------------------------------------------------------------------- %

% Matrix Size
n = 5;

% Entries for matrices
population = [-1, -1/10000, 0, 1/10000, 1];

% The generator
g = @() randomMatrix(population, n);

% Set the working directory
workingDir = '~/Real5x5_Error/';

% Set the options
margin = struct('bottom', -2.5, ...
                   'top',  2.5, ...
                  'left', -2.5, ...
                 'right',  2.5);

opts = struct(   'numDataFiles', 5,      ...
              'matricesPerFile', 1e6,    ...
                       'height', 1001,   ...
                       'margin', margin, ...
                   'ignoreReal', true);

% Generate the eigenvalue data
generateRandomSample(g, workingDir, opts);

% Process the data
fname = processData(workingDir, opts);

% The colormap
T = [ 34,    0,   0;
      78,    0,   0;
      139,   0,   0;
      226,   0,   0;
      236,  46,   0;
      246, 134,   0;
      254, 201,   0;
      255, 220,  57;
      255, 231, 119;
      255, 243, 185;
      255, 255, 255]/255;
x = [0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0];

% Make an image
processImage(workingDir, fname, T, x);
