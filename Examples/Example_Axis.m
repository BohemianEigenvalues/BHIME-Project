% ---------------------------------------------------------------------------- %
% This example will generate a plot in the matlab plot viewer window that      %

% contains axis labels.                                                        %
% ---------------------------------------------------------------------------- %
workingDir = '~/Desktop/Example_Axis/';

g = @() randomMatrix([-1, 0, 1], 5);

margin = struct('bottom', -3.3, ...
                   'top',  3.3, ...
                  'left', -5,   ...
                 'right',  5);

opts = struct('numDataFiles', 1, ...
           'matricesPerFile', 1e6, ...
                    'margin', margin, ...
                    'height', 2001, ...
                'symmetryIm', true, ...
            'backgroundColor', [1, 1, 1]);

generateRandomSample(g, workingDir, opts);

pFilename = processData(workingDir, opts);

% This will load a variable called "mesh" that is a matrix of integer counts of
% the eigenvalues at each point (pixel).
load(fullfile(workingDir, 'ProcessedData', pFilename));

% This is used by the image function to set the axes
xx = linspace(margin.left, margin.right, size(mesh, 2));
yy = linspace(margin.bottom, margin.top, size(mesh, 1));

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

% See processImage.m lines 159-188

% Find the maximum (non-infinite) density
maxDensity = max(mesh(isfinite(mesh)));

% Find the minimum (non-zero) density
minDensity = min(mesh(mesh ~= 0));

% Take the logarithm of the densities
mesh = log(double(mesh));

% Convert max/min densities to log scale
bounds = [log(double(minDensity)), log(double(maxDensity))];

% Interpolate the colormap at 255 values
map = interp1(x, T, linspace(0, 1, 255));
    
% Convert to rgb based on the color map
rgb = double2rgb(mesh, map, bounds, opts.backgroundColor);

% Use the image function to plot a color image
figure();
image(xx, yy, rgb);
set(gca,'YDir','normal');  % Flip the plot so it isn't upsidedown
axis equal;
