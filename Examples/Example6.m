% ----------------------------------------------------------------------- %
% EXAMPLE 6 -- Eigenvalues of 5x5 matrices with entries in                %
%              {exp(2*pi*1i*k/5) | k = 0...4}                             %
% ----------------------------------------------------------------------- %

% Matrix Size
n = 5;

% Entries for matrices
population = exp(2*pi*1i*(0:4)/5);

% The generator
g = @() randomMatrix(population, n);

% Set the working directory
workingDir = '~/UnitCircle_5/';

% Set the options
margin = struct('bottom', -3.5, ...
                   'top',  3.5, ...
                  'left', -3.5, ...
                 'right',  3.5);

map = @(z) [z;
            abs(z).*exp((angle(z) + 2*pi/5)*1i);
            abs(z).*exp((angle(z) + 4*pi/5)*1i);
            abs(z).*exp((angle(z) + 6*pi/5)*1i);
            abs(z).*exp((angle(z) + 8*pi/5)*1i)];
             
opts = struct(   'numDataFiles', 5,      ...
              'matricesPerFile', 1e7,    ...
                       'height', 2001,   ...
                       'margin', margin, ...
                          'map', map);

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
x = [0, 0.1, 0.18, 0.24, 0.28, 0.35, 0.4, 0.48, 0.55, 0.6, 1.0];

% Make an image
processImage(workingDir, fname, T, x);
