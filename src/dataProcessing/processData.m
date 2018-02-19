% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Feb. 19/2018                                                %
%                                                                         %
% This function will read all .mat files created by the                   %
% generateRandomSample function and create a histogram in the complex     %
% plane of their densities.                                               %
%                                                                         %
%   - A readme file will be automatically generated the first time this   %
%     function is called, if a readme file already exists then it will    %
%     be appended with information about this call.                       %
%   - A new directory ProcessedData will be created in the workingDir     %
%   - This function will create temporary files in the directory          %
%     containing the data while it is running.                            %
%                                                                         %
% INPUT                                                                   %
%   workingDirIn ... (str) The directory where files should be written    %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
%   filenamePrefix ............. Default = BHIME                          %
%                                The prefix for the .mat files that       %
%                                contain the eigenvalues.                 %
%   height ..................... Default = 1001 (pixels)                  %
%                                The height (in pixels) of the grid to be %
%                                used,  the width is determined from the  %
%                                margin such that each grid point is      %
%                                square.                                  %
%   ignoreReal ................. Default = false                          %
%                                Ignore any points that have zero         %
%                                imaginary part.                          %
%   ignoreRealTol .............. Default = 1e-8                           %
%                                The tolerance for how close the          %
%                                imaginary part of and eigenvalue is      %
%                                before it is considered to be a real     %
%                                number.                                  %
%   map ........................ Default = z -> z (no mapping)            %
%                                Map the eigenvalues by a given function  %
%                                handel. The function MUST be vectorized. %
%   margin ..................... Default = large enough to fit all the    %
%                                          points                         %
%                                A struct with keys:                      %
%                                    bottom, top, left, right             %
%                                that indicated the margins for the image %
%                                if more than one files is used it will   %
%                                only fit all the data in the data first  %
%                                file.                                    %
%   numProcessFiles ............ Default = All files in data directory    %
%                                The number of data files to process      %
%   outputFileType ............. Default = mat                            %
%                                Type of file to be output (either mat or %
%                                txt)                                     %
%   overrideDataDir ............ Default = [empty string] (not set)       %
%                                If the data is in a directory that is    %
%                                not workingDir/Data/ then set this to    %
%                                the directory containing the data files. %
%   overrideProcessedDataDir ... Default = [empty string] (not set)       %
%                                Set this option to a non-empty string    %
%                                indicating a directory to write the      %
%                                processed data files to. If not set they %
%                                will be written to the directory         %
%                                workingDir/ProcessedData/. If the        %
%                                directory does not exist the this        %
%                                function will create it.                 %
%   symmetryRe ................. Default = false                          %
%                                If true, symmetry across the real axis   %
%                                will be used to effectively double the   %
%                                number of points.                        %
%   symmetryIm ................. Default = false                          %
%                                If true, symmetry across the imaginary   %
%                                axis will be used to effectively double  %
%                                the number of points.                    %
%                                                                         %
% OUTPUT                                                                  %
%   fname ... A string of the name of file that data is written to        %
%   stats ... A struct containing statistics about the data:              %
%                 numUniquePts ... Number of points containing at least   %
%                                  one eigenvalue                         %
%                 outsideCount ... Number of points outside the margin    %
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
function [outputFilename, stats] = processData(workingDirIn, options)
    
    tic
    
    narginchk(1, 2)
    
    if nargin < 2
        options = struct();
    end
    
    % Process the options
    opts = processOptions(options);
    
    % Get the working directory
    if workingDirIn(end) ~= filesep
        workingDir = [workingDirIn, filesep];
    else
        workingDir = workingDirIn;
    end
    
    % Get the Data directory
    if opts.overrideDataDirIsSet
        dataDir = opts.overrideDataDir;
    else
        dataDir = [workingDir, 'Data', filesep];
    end
    
    % Get and create the ProcessedData directory
    if opts.overrideProcessedDataDirIsSet
        processedDataDir = opts.overrideProcessedDataDir;
    else
        processedDataDir = [workingDir, 'ProcessedData', filesep];
    end
    mkdir_if_not_exist(processedDataDir);
    
    % Get numProcessFiles
    if ~opts.numProcessFilesIsSet
        opts.numProcessFiles = getnumProcessFiles(dataDir, opts.filenamePrefix);
    end
    
    % Get default margin
    if ~opts.marginIsSet
        opts.margin = getDefaultMargin(dataDir, opts.filenamePrefix);
    end
    
    % Get the resolution
    opts.resolution = getResolution(opts.margin, opts.height);
    
    % Output filename
    outputFilename = makeOutputFilename(processedDataDir, opts);
    
    % Write readme file
    writeReadMe(processedDataDir, outputFilename, opts);
    
    [mesh, stats] = process_density(dataDir, opts);
    
    % Write mesh to a text file
    if strcmp(opts.outputFileType, 'txt')
        dlmwrite([processedDataDir, outputFilename], mesh, 'delimiter',' ','precision',15);
    elseif strcmp(opts.outputFileType, 'mat')
        save([processedDataDir, outputFilename], 'mesh', 'stats', '-v7.3');
    end
    
    toc
end


% ======================================================================= %
% FUNCTIONS                                                               %
% ======================================================================= %


% ----------------------------------------------------------------------- %
% process_density                                                         %
%                                                                         %
% Process all data files.                                                 %
%                                                                         %
% INPUT                                                                   %
%   dataDir ... Directory containing data files                           %
%   opts ...... Options struct                                            %
%                                                                         %
% OUTPUT                                                                  %
%   totalMesh ... The cumulative mesh for all files                       %
%   stats ....... Statistics about the data                               %
% ----------------------------------------------------------------------- %
function [totalMesh, stats] = process_density(dataDir, opts)
    
    % Get the options
    resolution      = opts.resolution;
    margin          = opts.margin;
    filenamePrefix  = opts.filenamePrefix;
    numProcessFiles = opts.numProcessFiles;
    
    % Make the mesh to store the result
    totalMesh = uint64(zeros(resolution.height, resolution.width));
    
    % Sizes of the points
    pointWidth  = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;
    
    % Process all the data files
    parfor i = 1:numProcessFiles
        
        localMesh = uint64(zeros(resolution.height, resolution.width));
        
        fprintf('File %d of %d\n', i, numProcessFiles);
        
        dataFilename = [dataDir, filenamePrefix, ...
                        '_', num2str(i), '.mat'];
        
        % Load the eigenvalues
        z = parLoad(dataFilename);
        
        % Check if there are weights for the eigenvalues
        hasWeights = ismember('weights', fieldnames(z));
        if hasWeights
            weights = z.weights(:);
        else
            weights = 1;    % Just to make Matlab happy :/
        end
        
        % Get the eigenvalues
        z = z.eigVals(:);
        
        % Map the eigenvalues
        z = opts.map(z);
        
        % If symmetry, add reflection of values across imaginary axis
        if opts.symmetryIm
            valid = abs(imag(z)) > eps;    % Ensure points aren't doubled
            z = [z; -conj(z(valid))];
            
            % Update weights
            if hasWeights
                weights = [weights; weights(valid)];
            end
            
        end
        if opts.symmetryRe
            valid = abs(real(z)) > eps;    % Ensure points aren't doubled
            z = [z; conj(z(valid))];
            
            % Update weights
            if hasWeights
                weights = [weights; weights(valid)];
            end
            
        end
        
        % Remove points not in margin
        valid = isInMargin(real(z), imag(z), margin);
        z = z(valid);
        
        % If there are valid points
        if length(z) ~= 0
            
            % Update weights
            if hasWeights
                weights = weights(valid);
            end
            
            % tolerance for ignoreReal option
            tol = opts.ignoreRealTol;
            
            % Ignore real points if ignoreReal option is true
            if opts.ignoreReal
                valid = abs(imag(z)) > tol
                z = z(valid);
                
                % Update weights
                if hasWeights
                    weights = weights(valid);
                end
            end
            
            xVal = uint32(ceil((real(z) - margin.left)/pointWidth));
            yVal = uint32(ceil((imag(z) - margin.bottom)/pointHeight));
            
            idx = uint32(((xVal - 1)*resolution.height + yVal));
            
            if hasWeights
                for j=1:length(weights)
                    localMesh(idx(j)) = localMesh(idx(j)) + uint64(weights(j));
                end
            else
                % Count the occurences of each unique value
                y = sort(idx);
                p = find([true; diff(y)~=0; true]);
                values = y(p(1:end-1));
                instances = diff(p);
                localMesh(values) = uint64(instances);
            end
            
            totalMesh = totalMesh + localMesh;
        end
        
    end
    
    stats = struct();
    
end


% ----------------------------------------------------------------------- %
% writeReadMe                                                             %
%                                                                         %
% Write a readme file in the processedDataDir with information about when %
% the data was created, what was used to create the data, etc.            %
% ----------------------------------------------------------------------- %
function writeReadMe(processedDataDir, outputFilename, opts)
    
    % Get the options
    ignoreReal      = opts.ignoreReal;
    ignoreRealTol   = opts.ignoreRealTol;
    map             = opts.map;
    margin          = opts.margin;
    numProcessFiles = opts.numProcessFiles;
    resolution      = opts.resolution;
    symmetryIm      = opts.symmetryIm;
    symmetryRe      = opts.symmetryRe;
    
    % Write a readme file
    file = fopen([processedDataDir, 'README.txt'],'a');
    
    % File name
    fprintf(file, [outputFilename, '\n']);
    
    % Date
    fprintf(file, '    Created ............. %s\n', ...
                   datestr(now,'mmmm dd/yyyy HH:MM:SS AM'));
    
    % ignoreReal
    fprintf(file, '    ignoreReal .......... ');
    if ignoreReal
        fprintf(file, 'true\n');
        fprintf(file, '    ignoreRealTol ....... %.5E\n', ...
                       ignoreRealTol);
    else
        fprintf(file, 'false\n');
    end
    
    % Map
    fprintf(file, '    map ................. %s\n', func2str(map));
    
    % Margin
    fprintf(file, '    margin .............. bottom: %.5f\n', ...
                   margin.bottom);
    fprintf(file, '                             top: %.5f\n', ...
                   margin.top);
    fprintf(file, '                            left: %.5f\n', ...
                   margin.left);
    fprintf(file, '                           right: %.5f\n', ...
                   margin.right);
    
    % numProcessFiles
    fprintf(file, '    numProcessFiles ..... %d\n', numProcessFiles);
    
    % Resolution
    fprintf(file, '    resolution .......... %d x %d\n', ...
                   resolution.width, resolution.height);
    
    % SymmetryIm
    fprintf(file, '    symmetryIm .......... ');
    if symmetryIm
        fprintf(file, 'true\n');
    else
        fprintf(file, 'false\n');
    end
    
    % SymmetryRe
    fprintf(file, '    symmetryRe .......... ');
    if symmetryRe
        fprintf(file, 'true\n');
    else
        fprintf(file, 'false\n');
    end
    
    fprintf(file, '\n\n\n');
    fclose(file);
    
end

% -------------------------------------------------------------------------
function b = isInMargin(x, y, margin)
    
    b = x > margin.left & x < margin.right & y > margin.bottom & y < margin.top;
    
end


% -------------------------------------------------------------------------
% Use in place of load in a parfor loop
function l = parLoad(fName)
    l = load(fName);
end


% -------------------------------------------------------------------------
% Use in place of save in a parfor loop
function parSave(fname, numvars, varargin)
    for i = 1:numvars
       eval([inputname(i+2),' = varargin{i};']);  
    end
    save('-mat',fname,inputname(3));
    for i = 2:numvars    
        save('-mat',fname,inputname(i+2),'-append');
    end
end


% ----------------------------------------------------------------------- %
% getnumProcessFiles                                                      %
%                                                                         %
% Determine the index of the last data file                               %
%                                                                         %
% OUTPUT                                                                  %
%   Integer, the index of thelast file in the data directory with         %
%   the filenamePrefix prefix.                                            %
% ----------------------------------------------------------------------- %
function n = getnumProcessFiles(dataDir, filenamePrefix)
    
    folderInfo = dir(dataDir);
    
    n = 0;
    
    dfpLen = numel(filenamePrefix);
    
    for k=1:length(folderInfo)
        
        filename = folderInfo(k).name;
        [~, name, ext] = fileparts(filename);
        
        if ~strcmp(ext, '.mat')
            continue
        end
        if numel(name) < dfpLen + 2
            continue
        end
        if ~strcmp(name(1:dfpLen+1), [filenamePrefix, '_'])
            continue
        end
        
        i = str2double(name(dfpLen+2:end));
        
        n = max(i, n);
    end
    
end


% ----------------------------------------------------------------------- %
% getDefaultMargin                                                        %
%                                                                         %
% Determine the default margin such that all data in first file fits      %
% in margin                                                               %
%                                                                         %
% OUTPUT                                                                  %
%   Margin struct                                                         %
% ----------------------------------------------------------------------- %
function m = getDefaultMargin(dataDir, dataFilePrefix)
    
    % Read the first file
    dataFilename = [dataDir, dataFilePrefix, '_1.mat'];
    data = parLoad(dataFilename);
    
    eigVals = data.eigVals;
    
    bottom = min(imag(eigVals(:)));
    top    = max(imag(eigVals(:)));
    left   = min(real(eigVals(:)));
    right  = max(real(eigVals(:)));
    
    m = struct('bottom', bottom, ...
               'top', top, ...
               'left', left, ...
               'right', right);
    
end


% ----------------------------------------------------------------------- %
% outputFilename                                                          %
%                                                                         %
% Determine the name of the output file.                                  %
%                                                                         %
% OUTPUT                                                                  %
%   A string of the form                                                  %
%   {dataFilePrefix}-{width}x{height}-{colorBy}-{symmetry}                %
% ----------------------------------------------------------------------- %
function outputFilename = makeOutputFilename(processedDataDir, opts)
    
    resolution      = opts.resolution;
    filenamePrefix  = opts.filenamePrefix;
    numProcessFiles = opts.numProcessFiles;
    outputFileType  = opts.outputFileType;
    
    outPrefix = [filenamePrefix, '-', ...
                 num2str(resolution.width), 'x', ...
                 num2str(resolution.height)];

    if exist([processedDataDir, outPrefix, '.', outputFileType]) == 2
        outPrefix = [outPrefix, '-'];

        i = 2;
        while exist([processedDataDir, outPrefix , num2str(i), '.', outputFileType]) == 2
            i = i + 1;
        end
        outputFilename = [outPrefix, num2str(i), '.'];
    else
        outputFilename = [outPrefix, '.'];
    end
    outputFilename = [outputFilename, outputFileType];
end


% ----------------------------------------------------------------------- %
% getResolution                                                           %
%                                                                         %
% Compute the resolution struct based on the height and margins.          %
%                                                                         %
% OUTPUT                                                                  %
%   A struct resolution = {'width', w, 'height', h}                       %
% ----------------------------------------------------------------------- %
function resolution = getResolution(margin, height)
    
    % Check the margins and make the resolution structure
    if margin.bottom >= margin.top
        error 'Bottom margin must be less than top margin';
    end
    if margin.left >= margin.right
        error 'Left margin must be less than top margin';
    end
    
    width = getWidth(margin, height);
    
    resolution = struct('width', width, 'height', height);
    
end


% ----------------------------------------------------------------------- %
% getWidth                                                                %
%                                                                         %
% Compute the width (in px) based on the height and the margins such      %
% that each grid point is a square                                        %
%                                                                         %
% OUTPUT                                                                  %
%   A struct resolution = {'width', w, 'height', h}                       %
% ----------------------------------------------------------------------- %
function width = getWidth(margin, height)
    heightI = margin.top - margin.bottom;
    widthI = margin.right - margin.left;
    width = floor(widthI*height/heightI);
end