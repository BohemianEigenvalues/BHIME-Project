% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2018)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Oct. 28/2018                                                %
%                                                                         %
% Convert processed data into an image given a colormap. A new directory  %
% Images will be created in the working directory if it doesn't already   %
% exist.                                                                  %
%                                                                         %
% INPUT                                                                   %
%   processedDataFilename ... (str) Name of the file containing the data  %
%                             to be plotted                               %
%   workingDir .............. (str) The directory where files should be   %
%                             written                                     %
%   T ....................... (optional) A m x 3 matrix where each row    %
%                             represents RGB values                       %
%   x ....................... (optional) Weights for the color map        %
%                             (vector of m values)                        %
%                                                                         %
% OPTIONS                                                                 %
%   backgroundColor ............ Vector with 3 elements specifying RGB    %
%                                values to use for the background color,  %
%                                values should be in [0,1] (i.e. [0,0,0]  %
%                                is black, [1,1,1] is white).             %
%   maxDensity ................. Default = Largest finite density point   %
%                                in the processed data file.              %
%                                Maximum density for coloring the image.  %
%                                If any points have a higher density,     %
%                                they will be set to this value.          %
%   maxDensity ................. Default = Largest non-zero density point %
%                                in the processed data file.              %
%                                Minimum density for coloring the image.  %
%                                If anypoints have a lower density        %
%                                (excluding zero), they will be set to    % 
%                                this value.                              %
%   overrideProcessedDataDir ... Default = [empty string] (not set)       %
%                                Set this option to a non-empty string    %
%                                indicating a directory to write the      %
%                                processed data files to. If not set they %
%                                will be written to the directory         %
%                                workingDir/ProcessedData/. If the        %
%                                directory does not exist the this        %
%                                function will create it.                 %
%   overrideImagesDir .......... Default = [empty string] (not set)       %
%                                Set this option to a non-empty string    %
%                                indicating a directory to write the      %
%                                images to. If not set they will be       %
%                                written to the directory                 %
%                                workingDir/Images/. If the directory     %
%                                does not exist this function will        %
%                                create it.                               %
%                                                                         %
% OUTPUT                                                                  %
%   outputImageFilename ... The name of the image file that will be       %
%                           written                                       %
%                                                                         %
% LICENSE                                                                 %
%   This program is free software: you can redistribute it and/or modify  %
%   it under the terms of the GNU General Public License as published by  %
%   the Free Software Foundation, either version 3 of the License, or     %
%   any later version.                                                    %
%                                                                         %
%   This program is distributed in the hope that it will be useful,       %
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%   GNU General Public License for more details.                          %
%                                                                         %
%   You should have received a copy of the GNU General Public License     %
%   along with this program.  If not, see http://www.gnu.org/licenses/.   %
% ----------------------------------------------------------------------- %
function outputImageFilename = processImage(processedDataFilename, workingDirIn, cmap, x, options)
    
    narginchk(2, 5);
    
    % processImage(processedDataFilename, workingDirIn)
    if nargin == 2
        options = struct();
        cmap = [  0,   0,   0;
                  0,   0, 255;
                  0, 255, 255;
                  0, 255,   0;
                255, 255,   0;
                255,   0,   0;
                255,   0,   0;
                255, 255, 255;
                255, 255, 255]/255;
        x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];
    end
    
    % processImage(processedDataFilename, workingDirIn, options)
    if nargin == 3
        cmap = [  0,   0,   0;
                  0,   0, 255;
                  0, 255, 255;
                  0, 255,   0;
                255, 255,   0;
                255,   0,   0;
                255,   0,   0;
                255, 255, 255;
                255, 255, 255]/255;
        x = [0, 0.1, 0.16, 0.22, 0.28, 0.34, 0.4, 0.55, 1.0];
    end
    
    % processImage(processedDataFilename, workingDirIn, cmap, x)
    if nargin == 4
        options = struct();
    end
    
    % Check that cmap and x have matching dimension
    if size(cmap, 2) ~= 3
        error('Color map must have 3 columns');
    end
    if size(cmap, 1) ~= max(size(x))
        error('Colormap and weight vector must be of same length');
    end
    
    % Check that the values in cmap are in [0, 1]
    if max(cmap(:)) > 1 || min(cmap(:)) < 0
        error('Color map must have values in [0, 1]');
    end
    
    % Process the options
    opts = processOptions(options);
    
    % Get the working directory
    if workingDirIn(end) ~= filesep
        workingDir = [workingDirIn, filesep];
    else
        workingDir = workingDirIn;
    end
    
    % Get the ProcessedData directory
    if opts.overrideProcessedDataDirIsSet
        processedDataDir = opts.overrideProcessedDataDir;
    else
        processedDataDir = [workingDir, 'ProcessedData', filesep];
    end
    
    % Get the Images directory
    if opts.overrideImagesDirIsSet
        imagesDir = opts.overrideImagesDir;
    else
        imagesDir = [workingDir, 'Images', filesep];
    end
    mkdir_if_not_exist(imagesDir);
    
    % The file to convert into an image
    dataFile = [processedDataDir, processedDataFilename];
    
    % Load the input file
    load(dataFile, 'mesh');
    
    % Make the name for the output file
    outputImageFilename = makeOutputImageFilename(imagesDir);
    fprintf('Image will be written to: %s\n', outputImageFilename);
    
    % ------------------------
    
    % Set the maximum value
    if opts.maxDensityIsSet
        maxDensity = opts.maxDensity;
        mesh = min(mesh, maxDensity);
    else
        maxDensity = max(mesh(isfinite(mesh)));
    end
    
    % Set the minimum value
    if opts.minDensityIsSet
        minDensity = opts.minDensity;
        mesh(mesh ~= 0) = max(mesh(mesh ~= 0), minDensity);
    else
        minDensity = min(mesh(mesh ~= 0));
    end
    
    % Take the logarithm of the densities
    mesh = log(double(mesh));
    
    % Convert max/min densities to log scale
    bounds = [log(double(minDensity)), log(double(maxDensity))];
    
    % Flip so image isn't upside down
    mesh = flipud(mesh);
    
    % Interpolate the colormap at 255 values (why?)
    map = interp1(x, cmap, linspace(0, 1, 255));
    
    % Convert to rgb based on the color map
    rgb = double2rgb(mesh, map, bounds, opts.backgroundColor);
    
    % Save the image
    imwrite(rgb, [imagesDir, outputImageFilename]);
    
    % Write a readme
    writeReadMe(imagesDir, ...
                outputImageFilename, ...
                processedDataFilename, ...
                minDensity, ...
                maxDensity, ...
                cmap, x, ...
                opts);
    
end


% ----------------------------------------------------------------------- %
% outputFilename                                                          %
%                                                                         %
% Determine the name of the output file.                                  %
%                                                                         %
% OUTPUT                                                                  %
%   A string of the form                                                  %
%   Image-k.png for an positive integer k                                 %
% ----------------------------------------------------------------------- %
function outputFilename = makeOutputImageFilename(imagesDir)
    
    outPrefix = ['Image-'];
    
    k = 1;
    while exist([imagesDir, outPrefix, num2str(k), '.png']) == 2
        k = k + 1;
    end
    outputFilename = [outPrefix, num2str(k), '.png'];
    
end


% ----------------------------------------------------------------------- %
% writeReadMe                                                             %
%                                                                         %
% Write a readme file in the processedDataDir with information about when %
% the data was created, what was used to create the data, etc.            %
% ----------------------------------------------------------------------- %
function writeReadMe(imagesDir, outputImageFilename, processedDataFilename, minDensity, maxDensity, cmap, x, opts)
    
    % Append to readme file
    file = fopen([imagesDir, 'README.txt'], 'a');
    
    % File name
    fprintf(file, [outputImageFilename, '\n']);
    
    % Date
    fprintf(file, '    created ................. %s\n', ...
                   datestr(now,'mmmm dd/yyyy HH:MM:SS AM'));
    
    % backgroundColor
    fprintf(file, '    backgroundColor ......... [%.5f, %.5f, %.5f]\n', opts.backgroundColor(1), opts.backgroundColor(2), opts.backgroundColor(3));
    
    % minDensity
    fprintf(file, '    minDensity .............. %u\n', uint64(minDensity));
    
    % maxDensity
    fprintf(file, '    maxDensity .............. %u\n', uint64(maxDensity));
    
    % processedDataFilename
    fprintf(file, '    processedDataFilename ... %s\n', processedDataFilename);
    
    % Colormap
    fprintf(file, '    colormap ................ ');
    fprintf(file, '[%.5f, %.5f, %.5f]\n', cmap(1, 1), cmap(1, 2), cmap(1, 3));
    for i=2:size(cmap, 1)
        fprintf(file, '                              ');
        fprintf(file, '[%.5f, %.5f, %.5f]\n', cmap(i, 1), cmap(i, 2), cmap(i, 3));
    end
    
    % Colormap weights
    fprintf(file, '    colormap weights ........ ');
    fprintf(file, '[');
    for i=1:length(x)-1
        fprintf(file, '%.5f, ', x(i));
    end
    fprintf(file, '%.5f, ', x(end));
    fprintf(file, ']\n\n');
    
    fclose(file);
    
end