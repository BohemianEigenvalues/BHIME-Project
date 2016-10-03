% ----------------------------------------------------------------------- %
% EXAMPLE 3 -- Eigenvalues of 5x5 matrices with entries in                %
%              {-20, -1, 0, 1, 20}                                        %
% ----------------------------------------------------------------------- %

% Matrix size
n = 5;

% Entries for matrices
population = [-20, -1, 0, 1, 20];

% The generator
g = @() randomMatrix(population, n);

workingDir = '~/Real5x5_wide/';

% Set the options
margin = struct('bottom', -45, ...
                   'top',  45, ...
                  'left', -55, ...
                 'right',  55);
            
% The ignoreReal option will ignore any values with zero complex part
% in the processData set. That is, values on the real line will be ignored.
% This can be useful for visualization as the density on the real axis is 
% typically much higher than elsewhere. The ignoreRealTol options
% specifies how close values must be to the real axis to be considered real
% values. i.e. abs(Im(z)) <  ignoreRealTol => z is real and hence ignored.
opts = struct('numDataFiles', 10, ...
           'matricesPerFile', 1e6, ...
                    'height', 501, ...
                    'margin', margin, ...
                'ignoreReal', true, ...
             'ignoreRealTol', 1e-10);

% Generate the data (may take a few minutes)
generateRandomSample(g, workingDir, opts);

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
x = [0, 0.1, 0.2, 0.3,   0.4,   0.5,  0.6, 0.8, 1.0];

% Make an image
processImage(workingDir, fname, T, x);