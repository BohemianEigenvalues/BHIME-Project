% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Mar. 8/2016                                                 %
%                                                                         %
% Convert processed data from a .mat or .txt file into an image given a   %
% color map. A new directory Images will be created in the working        %
% directory if it doesn't already exist.                                  %
%                                                                         %
% INPUT                                                                   %
%   workingDir .............. (str) The directory where files should be   %
%                             written                                     %
%   processedDataFilename ... (str) Name of the file containing the data  %
%                             to be plotted                               %
%   T ....................... A m x 3 matrix where each column represents %
%                             RGB values                                  %
%   x ....................... Weights for the color map (1xm matrix)      %
%                                                                         %
% TO DO                                                                   %
%   - Add default color maps                                              %
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
function processImage(workingDir, processedDataFilename, T, x)
    
    % Check the number of arguments
    if nargin < 4
        error('processData:notEnoughInputArguments', ...
              'function requires at least 4 input values');
    elseif nargin > 4
        error('processData:tooManyInputArguments', ...
        'function requires at most 4 input values');
    end
    
    % Make Images directory if it doesn't exist
    if workingDir(end) == filesep
        imagesDir = [workingDir, 'Images', filesep];
        processDataDir = [workingDir, 'ProcessedData', filesep];
    else
        imagesDir = [workingDir, filesep, 'Images', filesep];
        processDataDir = [workingDir, filesep, 'ProcessedData', filesep];
    end
    mkdir_if_not_exist(imagesDir)
    
    dataFile = [processDataDir, processedDataFilename];
    
    % check if input file is a .txt of .mat file and load data
    if strcmp(dataFile(end-2:end), 'txt')
        mesh = load(dataFile);
    elseif strcmp(dataFile(end-2:end), 'mat')
        load(dataFile, 'mesh');
    else
        error 'Invalid data file';
    end
    
    % Make the name for the output file
    outputImageFilename = makeOutputImageFilename()
    disp(['File will be written to: ', outputImageFilename]);
    
     % Write a readme file
    file = fopen([imagesDir, 'README.txt'], 'a');
    fprintf(file, ['Last Update: ', datestr(now, 'mmmm dd/yyyy HH:MM:SS AM'), '\n\n']);
    fprintf(file, [outputImageFilename, '\n']);
    fprintf(file, ['Data from: ', processedDataFilename]);
    fprintf(file, '\n\n\n');
    fclose(file);
    
    % ------------------------
    % Interpolate data
    mesh = log(double(mesh));
    
    map = interp1(x', T'/255, linspace(0, 1, 255));
    % -----------------
    
    % Apply the color map
    rgb = double2rgb(mesh, map);
    
    % Write the image to a file
    imwrite(rgb, [imagesDir, outputImageFilename]);
    
    
    % ------------------------------------------------------------------- %
    % outputFilename                                                      %
    %                                                                     %
    % Determine the name of the output file.                              %
    %                                                                     %
    % OUTPUT                                                              %
    %   A string of the form                                              %
    %   Image-k.png for an positive integer k                             %
    % ------------------------------------------------------------------- %
    function outputFilename = makeOutputImageFilename()
        
        outPrefix = ['Image-'];
        
        k = 1;
        while exist([imagesDir, outPrefix, num2str(k), '.png']) == 2
            k = k + 1;
        end
        outputFilename = [outPrefix, num2str(k), '.png'];
        
    end
    
    
end