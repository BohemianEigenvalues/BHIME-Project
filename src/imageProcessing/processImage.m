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
% OPTIONS                                                                 %
%   backgroundColor ... Vector with 3 elements specifying RGB values to   %
%                       use for the background color (base 255)           %
%                                                                         %
% TO DO                                                                   %
%   - Make work with built in color maps without jet'*255                 %
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
function processImage(workingDir, processedDataFilename, cmap, x, options)
    
    % Check the number of arguments
    if nargin < 4
        error('processData:notEnoughInputArguments', ...
              'function requires at least 4 input values');
    elseif nargin > 5
        error('processData:tooManyInputArguments', ...
        'function requires at most 5 input values');
    elseif nargin == 4
        options = struct();
    end
    
    % Make Images directory if it doesn't exist
    if workingDir(end) == filesep
        imagesDir = [workingDir, 'Images', filesep];
        processDataDir = [workingDir, 'ProcessedData', filesep];
    else
        imagesDir = [workingDir, filesep, 'Images', filesep];
        processDataDir = [workingDir, filesep, 'ProcessedData', filesep];
    end
    mkdir_if_not_exist(imagesDir);
    
    % Process the options
    opts = processOptions();
    backgroundColor = opts.backgroundColor;
    
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
    outputImageFilename = makeOutputImageFilename();
    disp(['File will be written to: ', outputImageFilename]);
    
     % Write a readme file
    file = fopen([imagesDir, 'README.txt'], 'a');
    fprintf(file, ['Last Update: ', datestr(now, 'mmmm dd/yyyy HH:MM:SS AM'), '\n\n']);
    fprintf(file, [outputImageFilename, '\n']);
    fprintf(file, ['Data from: ', processedDataFilename]);
    fprintf(file, '\n\n\n');
    fclose(file);
    
    % ------------------------
    
    x = autoColormapWeights();
    
    % Interpolate data
    mesh = log(double(mesh));
    
    map = interp1(x, cmap, linspace(0, 1, 255));
    % -----------------
    
    % Apply the color map
    valid = isfinite(mesh); % Clean this up
    bounds = [min(mesh(valid)) max(mesh(valid))];   % Clean this up
    rgb = double2rgb(mesh, map, bounds, backgroundColor/255);
    
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
    
    
    % ------------------------------------------------------------------- %
    % autoColormapWeights                                                 %
    %                                                                     %
    % Determine the weights for the color map based on the distribution   %
    % of the data                                                         %
    %                                                                     %
    % OUTPUT                                                              %
    %   A vector of weights (from 0 to 1)                                 %
    % ------------------------------------------------------------------- %
    function X = autoColormapWeights()
        %{
        numColors = size(cmap,1);
        
        lmesh = mesh(isfinite(mesh));
        
        umesh = log(double(unique(lmesh)));
        umesh = umesh(isfinite(umesh));
        %semilogy(umesh);
        
        X = zeros(numColors, 1);
        
        X(1) = 0;
        
        for i=2:numColors
            
            idx = round(((i-1)/(numColors-1))*(length(umesh)-1));
            X(i) = umesh(idx+1);
            
        end
        
        X = X/X(end);   % Normalize
        %}
        % ------------------
        
        numColors = size(cmap, 1);
        
        lmesh = mesh(isfinite(mesh));
        
        lmesh_d = double(lmesh);
        
        [a, ~] = hist(lmesh_d, unique(lmesh_d));
        
        % YES!!!!!!
        %disp([a', unique(lmesh)]);
        
        %disp(log(a'));
        
        ap = a';
        ap = ap(2:end);
        
        ulmesh = unique(lmesh);
        %disp([ap, ulmesh(2:end), log(ap)]);
        
        %for i=2:length(ap)
        %    ap(i) = ap(i-1)+ap(i);
        %end
        
        lap = log(ap);
        %lap = ap;
        
        % Normalize lap to [0,1]
        lapmin = min(lap);
        lapmax = max(lap);
        disp(lapmin);
        disp(lapmax);
        
        lapN = (lap-lapmin);
        lapN = lapN/max(lapN);
        
        for i=2:length(lapN)
            lapN(i) = lapN(i-1)+lapN(i);
        end
        
        lapN = (lapN-min(lapN));
        lapN = lapN/max(lapN);
        
        disp(lapN);
        
        % Build X and get values from lapN
        X = zeros(numColors, 1);
        
        for i=1:numColors
            
            idx = round(((i-1)/(numColors-1))*(length(lapN)-1))+1;
            X(i) = lapN(idx);
            
        end
        X(1) = 0;
        X(end) = 1;
        
        % Make X unique...
        
        % End YES!!!!!!!
        % -------------------------
        
    end
    
    
    % ------------------------------------------------------------------- %
    % processOptions                                                      %
    %                                                                     %
    % Process the options input options struct. If an option is not in    %
    % the options struct the default value is used.                       %
    %                                                                     %
    % INPUT                                                               %
    %   options ... (struct) contains keys corresponding to the options   %
    %                                                                     %
    % OUTPUT                                                              %
    %   A struct opts with keys                                           %
    %       backgroundColor                                               %
    % TO DO                                                               %
    %   - Add type checking for options                                   %
    % ------------------------------------------------------------------- %
    function opts = processOptions()
        
        % Check that options is a struct
        if ~isstruct(options)
            error('processData:InvalidOptionsStruct', ...
                  'options argument must be a structured array');
        end
        
        fnames = fieldnames(options);
        
        optionNames = struct('backgroundColor', 'backgroundColor');
        
        if ~all(ismember(fnames, fieldnames(optionNames)))
            error('processData:InvalidOption',  ...
                  'Invalid option provided');
        end
        
        % Default values
        backgroundColor = [0, 0, 0];
        
        
        % height
        if isfield(options, 'backgroundColor')
            backgroundColor = options.backgroundColor;
        end
        
        opts = struct('backgroundColor', backgroundColor);
        
    end
    
end