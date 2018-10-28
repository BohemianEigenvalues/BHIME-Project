% ----------------------------------------------------------------------- %
% EXAMPLE 7 -- Eigenvalues of 5x5 matrices where the entries are sampled  %
%              from the set of Fibonacci numbers up to 514,229.           %
% ----------------------------------------------------------------------- %

workingDir = '~/Fibonacci/';

margin = struct('bottom', -550000, ...
                   'top',  550000, ...
                  'left', -550000,   ...
                 'right',  550000);

opts = struct('numDataFiles', 50, ...
            'startFileIndex', 1, ...
           'matricesPerFile', 1e6, ...
                    'height', 2001, ...
                    'margin', margin, ...
                'symmetryIm', true, ...
                'ignoreReal', true);

matrixSize = 5;


% The generator
g = @() randomMatrix([0,1,2,3,5,8,13,21,34,55,89,144,233, 377, 610, 987, 1597, 2584, 4181,6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229], matrixSize);

% Genderate the random sample
generateRandomSample(g, workingDir, opts);

% Process the data
pFilename = processData(workingDir,opts);

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

x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% Make the image
processImage(pFilename, workingDir, T, x, opts);