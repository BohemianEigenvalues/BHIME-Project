% ----------------------------------------------------------------------- %
% EXAMPLE 2 -- Eigenvalues of 5x5 matrices with entries in {-1, 0, 1}     %
%                                                                         %
% This example is the same as Example 1, except options are used to       %
% override some of the default values. See Example 1 for a more detailed  %
% explaination of what each line of code is doing.                        %
% ----------------------------------------------------------------------- %

% Set the working directory
workingDir = '~/Real5x5_d3/';

% The generator (5x5 matrices with entries sampled from {-1, 0, 1})
g = @() randomMatrix([-1, 0, 1], 5);

% The margin sets the bounds for the image, it is a struct with the keys
% bottom, top, left and right
margin = struct('bottom', -3.3, ...
                   'top',  3.3, ...
                  'left', -5,   ...
                 'right',  5);

% The options are always stored in a struct where the keys are the option
% names and the values are the values for the option. Some options I use
% often are illustrated here. To see the full list of options as well as
% their descriptions see the README.
opts = struct('numDataFiles', 5, ...
           'matricesPerFile', 1e6, ...
                    'margin', margin, ...
                    'height', 2001, ...
                'symmetryIm', true, ...
           'backgroundColor', [1, 1, 1]);

% A random sample of 5 million (numDataFiles*matricesPerFile) matrices will  
% be generated and stored in 5 data files (each file will contain the 
% eigenvalues of 1 million matrices).
generateRandomSample(g, workingDir, opts);

% The grid for the densities will have a height of 2001 (pixels) and the 
% bounds for the grid are determined by the margin. The symmetryIm option
% will be used to reflect all the eigenvalues from the data file across
% the imaginary axis, effectively doubling the number of points plotted.
%
% Mathematical Note: The eigenvalues could be reflected across the real
%                    axis by using the symemtryRe option. But, a 
%                    consequence of the complex conjugate root theorem is 
%                    that these eigenvalues will already be in the image. 
%                    Therefore this symmetry does not add anything to the 
%                    image.
%                    Why does this happen? Beacuse the matrices are integer
%                    matrices, their characteristic polynomials will be 
%                    polynomials with integer coefficients. Thus, if a+bi 
%                    is a root of the polynomial (and hence an eigenvalue 
%                    of the original matrix) then by the complex conjugate 
%                    root theorem, a-bi will also be a root (and hence an 
%                    eigenvalue).
pFilename = processData(workingDir, opts);

% The colormap
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
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% The background color of the output image will be white ([1,1,1] = white)
processImage(pFilename, workingDir, T, x, opts);