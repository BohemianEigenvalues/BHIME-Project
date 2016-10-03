% ----------------------------------------------------------------------- %
% EXAMPLE 1 -- Eigenvalues of 5x5 matrices with entries in {-1, 0, 1}     %
%                                                                         %
% A simple introductory example for exploring the eigenvalues of a random %
% sample of 5x5 matrices where the entries are sampled uniformly from     %
% {-1, 0, 1}. An extensive description of how each method works is given  %
% here. No options are supplied to any of the methods. See the other      %
% examples for how to override the default options.                       %
% ----------------------------------------------------------------------- %

% Set the working directory, 3 folders will be created in this directory:
%   Data, ProcessedData and Images
workingDir = '~/Real5x5_d3/';

% The generator (5x5 matrices with entries sampled from {-1, 0, 1})
g = @() randomMatrix([-1, 0, 1], 5);

% The generateRandomSample function is used for generating a random sample
% of matrices, computing their eigenvalues and storing them. It will create
% the Data directory within the workingDir and write a .mat file  
% containing the eigenvalues it computes.
% The default size of the random sample is floor(1e6/matrixSize), in this
% case a sample of 200,000 matrices will be generated.
generateRandomSample(g, workingDir);

% The processData function will create the processedData directory within 
% the workingDir. It will process the raw eigenvalue data that was 
% generated by the generateRandomSample function. It will compute the 
% density of eigenvalues at each point in a grid over the complex plane 
% that will later become the output image.
% By default, the size of the grid (image) is 1001 (pixels) and the bounds 
% of the image are determined such that all eigenvalues (in the first data 
% file) fall inside the bounds.
pFilename = processData(workingDir);

% The colormap for the output image (should always be a m x 3 matrix of 
% values in [0, 1]
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

% The weights for the colormap, should be a vector of lenth m with strictly 
% increasing values from 0 to 1.
% Note: the densities in the output image are on a logarithmic scale.
x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];

% The processImage function will create the Images directory within the 
% workingDir and save an image there. The name of the image is generated
% automatically such that images are never overwritten.
processImage(workingDir, pFilename, T, x);