% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Apr. 22/2016                                                %
%                                                                         %
% This function will read all .mat files created by the                   %
% generateRandomSample function and sort the eigenvalues onto the complex %
% plane. If the colorBy input argument is set to Density, then a grid     %
% in the complex plane will count the number of eigenvalues that fall     %
% into each grid point. If Cond is specified then the average condition   %
% number of each eigenvalue that falls in each gridpoint will be          %
% computed.                                                               %
%   - A readme file will be automatically generated the first time this   %
%     function is called, if a readme file already exists then it will    %
%     be appended with information about this call.                       %
%   - A new directory ProcessedData will be created in the workingDir     %
%   - This function will create temporary files in the directory          %
%     containing the data while it is running.                            %
%                                                                         %
% INPUT                                                                   %
%   workingDir ... (str) The directory where files should be written      %
%   colorBy ...... Either 'density' or 'cond'                             %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
%   height ........... Default = 1001 (pixels)                            %
%                      The height (in pixels) of the grid to be used, the %
%                      width is determined from the margin such that each %
%                      grid point is square                               %
%   margin ........... Default = large enough to fit all the points       %
%                      A struct with keys:                                %
%                           bottom, top, left, right                      %
%                      that indicated the margins for the image           %
%                      if more than one files is used it will only fit    %
%                      all the data in the first file                     %
%   dataFilePrefix ... Default = BHIME                                    %
%                      The prefix for the .mat files that contain the     %
%                      eigenvalues and their condition numbers            %
%   outputFileType ... Default = mat                                      %
%                      Type of file to be output (either mat or txt)      %
%   symmetry ......... Default = false                                    %
%                      If true symmetry across the real and imaginary     %
%                      axes is used to give effectively quadrupling       %
%                      the number of points                               %
%   numFiles ......... Default = All files                                %
%                      The number of files to process                     %
%   map .............. Default = z -> z (no mapping)                      %
%                      Map the eigenvalues by a given function handel     %
%                      The function MUST be vectorized                    %
%                                                                         %
% OUTPUT                                                                  %
%   fname ... A string of the name of file that data is written to        %
%   stats ... A struct containing statistics about the data:              %
%                 numUniquePts ... Number of points containing at least   %
%                                  one eigenvalue                         %
%                 outsideCount ... Number of points outside the margin    %
%                                                                         %
% TO DO                                                                   %
%   - Separate symmetry into symmetry across real axis and symmetry       %
%     across imaginary axis                                               %
%   - Add option to ignore real axis                                      %
%   - Add option to map eigenvalues to a function (i.e. x -> 1/x)         %
%   - Can the symmetry and non-symmetric function be combined to make the %
%     program more compact?                                               %
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
function [outputFilename, stats] = processData(workingDir, colorBy, options)
    
    tic
    
    narginchk(2, 3)
    
    if nargin < 3
        options = struct();
    end
    
    % Make ProcessData directory if it doesn't exist
    if workingDir(end) == filesep
        processDataDir = [workingDir, 'ProcessedData', filesep];
        dataDir = [workingDir, 'Data', filesep];
    else
        processDataDir = [workingDir, filesep, 'ProcessedData', filesep];
        dataDir = [workingDir, filesep, 'Data', filesep];
    end
    mkdir_if_not_exist(processDataDir);
    
    % Process the options
    opts = processOptions(options);
    margin         = opts.margin;
    dataFilePrefix = opts.filenamePrefix;
    outputFileType = opts.outputFileType;
    symmetry       = opts.symmetry;
    numFiles       = opts.numProcessFiles;
    height         = opts.height;
    map            = opts.map;
    
    if ~opts.numProcessFilesIsSet
        numFiles = getNumFiles();
    end
    if ~opts.marginIsSet
        margin = getDefaultMargin()
    end
    
    % Get resolution
    resolution = getResolution();
    
    % Output filename
    outputFilename = makeOutputFilename();
    
    % Write readme file
    writeReadMe();
    
    % -------------------------
    
    % Call correct function for colorBy argument
    if strcmpi(colorBy, 'density')
        % Call the method for density
        [mesh, totalOutsideCount] = process_density();
    elseif strcmpi(colorBy, 'cond')
        % Call the method for condition number coloring
        [mesh, totalOutsideCount] = process_cond();
    else
        error 'Invalid value supplied for colorBy argument'
    end
    
    % -------------------------
    % Fill in the statistics struct
    stats = struct();
    
    % Number of unique points
    stats.numUniquePts = length(mesh(mesh ~= 0));
    stats.outsideCount = totalOutsideCount;
    
    % -------------------------
    
    % Write mesh to a text file
    if strcmp(outputFileType, 'txt')
        dlmwrite([processDataDir, outputFilename], mesh, 'delimiter',' ','precision',15);
    elseif strcmp(outputFileType, 'mat')
        save([processDataDir, outputFilename], 'mesh');
    end
    
    toc
    
    
    % =====================================================================
    % FUNCTIONS
    % =====================================================================
    
    
    % ---------------------------------------------------------------------
    function [mesh, totalOutsideCount] = process_density()
        if symmetry
            [mesh, totalOutsideCount] = process_density_symmetry();
        else
            [mesh, totalOutsideCount] = process_density_no_symmetry();
        end
    end
        
    % ---------------------------------------------------------------------
    
    function [mesh, totalOutsideCount] = process_cond()
        if symmetry
            [mesh, totalOutsideCount] = process_cond_symmetry();
        else
            [mesh, totalOutsideCount] = process_cond_no_symmetry();
        end
    end
    
    % ---------------------------------------------------------------------
    function [mesh, totalOutsideCount] = process_density_no_symmetry()
        
        totalOutsideCount = zeros(1, numFiles);
        
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename  = [processDataDir, 'tmp_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat'];
    
            % Call function to save the processed data to a temporary file
            totalOutsideCount(i) = process_density_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename, map);
            toc
    
        end
        
        totalOutsideCount = sum(totalOutsideCount);
        
        % ------------------------
    
        % Make the mesh to store the result (local to loop)
        mesh = uint32(zeros(resolution.height, resolution.width));
    
        % Read all temprary files and combine
        for i = 1:numFiles
    
            f = load([processDataDir,'tmp_',num2str(i),'.mat'],'mesh');
    
            mesh = mesh + f.mesh;
    
            delete([processDataDir,'tmp_',num2str(i),'.mat']);
    
        end
    
    end

    
    % ---------------------------------------------------------------------
    function [mesh, totalOutsideCount] = process_density_symmetry()
    
        totalOutsideCount = zeros(1, numFiles);
    
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename  = [processDataDir, 'tmp_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat'];
    
            % Call function to save the processed data to a temporary file
            totalOutsideCount = process_density_symmetry_tmp(resolution, margin, dataFilename, tmpFilename, map);
    
            toc
    
        end
        
        totalOutsideCount = sum(totalOutsideCount);
    
        % ------------------------
    
        % Make the mesh to store the result (local to loop)
        mesh = uint32(zeros(resolution.height, resolution.width));
    
        % Read all temprary files and combine
        for i = 1:numFiles
    
            f = load([processDataDir, 'tmp_',num2str(i),'.mat'],'mesh');
    
            mesh = mesh + f.mesh;
    
            delete([processDataDir, 'tmp_',num2str(i),'.mat']);
    
        end
    
    end
    
    
    % =====================================================================
    % =====================================================================
    
    function [mesh, totalOutsideCount] = process_cond_no_symmetry()
        
        totalOutsideCount = zeros(1, numFiles);
        
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename_count  = [processDataDir, 'tmp_count_',num2str(i),'.mat'];
            tmpFilename_total  = [processDataDir, 'tmp_total_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat']
    
            % Call function to save the processed data to a temporary file
            totalOutsideCount(i) = process_cond_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total, map);
    
            toc
    
        end
        
        totalOutsideCount = sum(totalOutsideCount);
        
        % ------------------------
    
        % Make the mesh to store the result
        count = uint32(zeros(resolution.height, resolution.width));
        total = zeros(resolution.height, resolution.width);
    
        % Read all temprary files and combine
        for i = 1:numFiles
    
            c = load([processDataDir,'tmp_count_',num2str(i),'.mat'],'count');
            t = load([processDataDir,'tmp_total_',num2str(i),'.mat'],'total');
    
            count = count + c.count;
            total = total + t.total;
    
            delete([processDataDir,'tmp_count_',num2str(i),'.mat']);
            delete([processDataDir,'tmp_total_',num2str(i),'.mat']);
    
        end
    
        % Compute the average
        mesh = zeros(resolution.height, resolution.width);
        
        for ii = 1:resolution.height
            for jj = 1:resolution.width
                if count(ii,jj) ~= 0
                    mesh(ii, jj) = total(ii, jj)/double(count(ii, jj));
    
                    if mesh(ii, jj) == 0 && total(ii, jj) ~= 0
                        disp(['underflow at (', num2str(ii), ', ', num2str(jj), ')']);
                        mesh(ii, jj) = realmin;
                    end
    
                end
            end
        end
    
    end
    
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    function [mesh, totalOutsideCount] = process_cond_symmetry()
    
        
        totalOutsideCount = zeros(1, numFiles);
        
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename_count  = [processDataDir, 'tmp_count_',num2str(i),'.mat'];
            tmpFilename_total  = [processDataDir, 'tmp_total_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat'];
    
            % Call function to save the processed data to a temporary file
            outsideCount(i) = process_cond_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total, map);
    
            toc
    
        end
        
        totalOutsideCount = sum(totalOutsideCount);
        
        % ------------------------
    
        % Make the mesh to store the result
        count = uint32(zeros(resolution.height, resolution.width));
        total = zeros(resolution.height, resolution.width);
    
        % Read all temprary files and combine
        for i = 1:numFiles
    
            c = load([processDataDir, 'tmp_count_',num2str(i),'.mat'],'count');
            t = load([processDataDir, 'tmp_total_',num2str(i),'.mat'],'total');
    
            count = count + c.count;
            total = total + t.total;
    
            delete([processDataDir, 'tmp_count_',num2str(i),'.mat']);
            delete([processDataDir, 'tmp_total_',num2str(i),'.mat']);
    
        end
    
        % Compute the average
        mesh = zeros(resolution.height, resolution.width);
        
        for ii = 1:resolution.height
            for jj = 1:resolution.width
                if count(ii,jj) ~= 0
                    mesh(ii, jj) = total(ii, jj)/double(count(ii, jj));
    
                    if mesh(ii, jj) == 0 && total(ii, jj) ~= 0
                        disp(['underflow at (', num2str(ii), ', ', num2str(jj), ')']);
                        mesh(ii, jj) = realmin;
                    end
    
                end
            end
        end
    
    end
    
    % ---------------------------------------------------------------------   
    
    % =====================================================================
    % MORE FUNCTIONS
    % =====================================================================
    
    
    % ------------------------------------------------------------------- %
    % getDefaultMargin                                                    %
    %                                                                     %
    % Determine the default margin such that all data in first file fits  %
    % in margin                                                           %
    %                                                                     %
    % OUTPUT                                                              %
    %   Margin struct                                                     %
    % ------------------------------------------------------------------- %
    function m = getDefaultMargin()
        
        % Read the first file
        dataFilename = [workingDir, 'Data', filesep, dataFilePrefix, '_1.mat'];
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
    
    
    % ------------------------------------------------------------------- %
    % getNumFiles                                                         %
    %                                                                     %
    % Determine the index of the last data file                           %
    %                                                                     %
    % OUTPUT                                                              %
    %   Integer, the index of thelast file in the data directory with     %
    %   the filenamePrefix prefix.                                        %
    % ------------------------------------------------------------------- %
    function n = getNumFiles()
        
        % Data directory
        if workingDir(end) == filesep
            dataDir = [workingDir, 'Data', filesep];
        else
            dataDir = [workingDir, filesep, 'Data', filesep];
        end
        
        folderInfo = dir(dataDir);
        
        n = 0;
        
        dfpLen = numel(dataFilePrefix);
        
        for k=1:length(folderInfo)
            
            filename = folderInfo(k).name;
            [~, name, ext] = fileparts(filename);
            
            if ~strcmp(ext, '.mat')
                continue
            end
            if numel(name) < dfpLen + 2
                continue
            end
            if ~strcmp(name(1:dfpLen+1), [dataFilePrefix, '_'])
                continue
            end
            
            i = str2double(name(dfpLen+2:end));
            
            n = max(i, n);
        end
        
    end
    
    
    % ------------------------------------------------------------------- %
    % getResolution                                                       %
    %                                                                     %
    % Compute the resolution struct based on the height and margins.      %
    %                                                                     %
    % OUTPUT                                                              %
    %   A struct resolution = {'width', w, 'height', h}                   %
    % ------------------------------------------------------------------- %
    function resolution = getResolution()
        
        % Check the margins and make the resolution structure
        if margin.bottom >= margin.top
            error 'Bottom margin must be less than top margin';
        end
        if margin.left >= margin.right
            error 'Left margin must be less than top margin';
        end
        
        width = getWidth();
        
        resolution = struct('width', width, 'height', height);
        
    end
    
    
    % ------------------------------------------------------------------- %
    % getWidth                                                            %
    %                                                                     %
    % Compute the width (in px) based on the height and the margins such  %
    % that each grid point is a square                                    %
    %                                                                     %
    % OUTPUT                                                              %
    %   A struct resolution = {'width', w, 'height', h}                   %
    % ------------------------------------------------------------------- %
    function width = getWidth()
        heightI = margin.top - margin.bottom;
        widthI = margin.right - margin.left;
        width = floor(widthI*height/heightI);
    end
    
    
    % ------------------------------------------------------------------- %
    % outputFilename                                                      %
    %                                                                     %
    % Determine the name of the output file.                              %
    %                                                                     %
    % OUTPUT                                                              %
    %   A string of the form                                              %
    %   {dataFilePrefix}-{width}x{height}-{colorBy}-{symmetry}            %
    % ------------------------------------------------------------------- %
    function outputFilename = makeOutputFilename()
        
        outPrefix = [dataFilePrefix, '-', num2str(resolution.width), 'x', num2str(resolution.height), '-', regexprep(colorBy,'(\<[a-z])','${upper($1)}')];
    
        if symmetry
            outPrefix = [outPrefix, '-sym-', num2str(numFiles)];
        else
            outPrefix = [outPrefix, '-', num2str(numFiles)];
        end
    
    
        if exist([processDataDir, outPrefix, '.', outputFileType]) == 2
            outPrefix = [outPrefix, '-'];
    
            i = 2;
            while exist([processDataDir, outPrefix , num2str(i), '.', outputFileType]) == 2
                i = i + 1;
            end
            outputFilename = [outPrefix, num2str(i), '.'];
        else
            outputFilename = [outPrefix, '.'];
        end
        outputFilename = [outputFilename, outputFileType];
    end
    
    
    % ------------------------------------------------------------------- %
    % writeReadMe                                                         %
    %                                                                     %
    % Write a readme file in the processDataDir with information about    %
    % when the data was created, what was used to create the data, etc.   %
    % ------------------------------------------------------------------- %
    function writeReadMe()
        
        file = fopen([processDataDir, 'README.txt'], 'a');
        
        fprintf(file, ['Last Update: ',datestr(now,'mmmm dd/yyyy HH:MM:SS AM'),'\n\n']);
        fprintf(file, [outputFilename, '\n']);
        fprintf(file, ['\tViewed on [', num2str(margin.left)]);
        if margin.bottom >= 0
            fprintf(file, '+');
        end
        fprintf(file, [num2str(margin.bottom), 'i, ', num2str(margin.right)]);
        if margin.top >= 0
            fprintf(file, '+');
        end
        fprintf(file, num2str(margin.top));
        fprintf(file, 'i]\n\n');
        
        if strcmpi(colorBy,'density')
            fprintf(file, '\tPlotted by density\n\n');
        elseif strcmpi(colorBy, 'cond')
            fprintf(file, '\tPlotted by eigenvalue condition number\n\n');
        end
        
        fprintf(file, ['\t', num2str(resolution.width), 'x', num2str(resolution.height), ' pixel grid']);
        
        fprintf(file, '\n\n\n');
        fclose(file);
        
    end
end





% -------------------------------------------------------------------------
% process_density_no_symmetry_tmp
% -------------------------------------------------------------------------
function outsideCount = process_density_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename, map)

    % Sizes of the points
    pointWidth = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;

    % Make the mesh to store the result (local to loop)
    mesh = uint32(zeros(resolution.height, resolution.width));
    
    % Number of points outside view
    outsideCount = uint64(0);

    % Load the eigenvalues
    eigVals = parLoad(dataFilename);
    eigVals = eigVals.eigVals;
    
    % Map the eigenvalues
    eigVals = map(eigVals);
    
    for i = 1:numel(eigVals)
        
        % Get the ith eigenvalue
        z = eigVals(i);
        
        % Get the real and imaginary parts
        xVal = real(z);
        yVal = imag(z);
        
        % 1 =====
        if isInMargin(xVal, yVal, margin)
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
        else
            outsideCount = outsideCount + 1;
        end
        
    end

    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(eigVals)*100);

    % Write mesh to a temporary .mat file
    parSave(tmpFilename, 1, mesh);

end


% -------------------------------------------------------------------------
function outsideCount = process_density_symmetry_tmp(resolution, margin, dataFilename, tmpFilename, map)

    % Sizes of the points
    pointWidth = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;

    % Make the mesh to store the result (local to loop)
    mesh = uint32(zeros(resolution.height, resolution.width));

    % Number of points outside view
    outsideCount = uint64(0);

    % Load the eigenvalues
    eigVals = parLoad(dataFilename);
    eigVals = eigVals.eigVals;
    
    % Map the eigenvalues
    eigVals = map(eigVals);

    for i = 1:numel(eigVals)

        % Get the ith eigenvalue
        z = eigVals(i);
        
        % Get the real and imaginary parts
        xVal = real(z);
        yVal = imag(z);
        
        % 1 =====
        if isInMargin(xVal, yVal, margin)
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
        else
            outsideCount = outsideCount + 1;
        end
        
        % 2 =====
        xVal = xVal * -1;
        if isInMargin(xVal, yVal, margin)
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
        else
            outsideCount = outsideCount + 1;
        end
        
        % 3 =====
        yVal = yVal * -1;
        if isInMargin(xVal, yVal, margin)
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
        else
            outsideCount = outsideCount + 1;
        end
        
        % 4 =====
        xVal = xVal * -1;
        if isInMargin(xVal, yVal, margin)
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
        else
            outsideCount = outsideCount + 1;
        end
    end

    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(eigVals)*100);

    % Write mesh to a temporary .mat file
    parSave(tmpFilename, 1, mesh);

end


% -------------------------------------------------------------------------
function outsideCount = process_cond_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total, map)
    
    % Sizes of the points
    pointWidth = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;
    
    % Make the mesh to store the result
    count = uint32(zeros(resolution.height, resolution.width));
    total = zeros(resolution.height, resolution.width);
    
    
    % Number of points outside view
    outsideCount = uint64(0);
    
    % Load the files
    data = parLoad(dataFilename);
    
    % Get the eigenvalues
    eigVals = data.eigVals;
    
    % Map the eigenvalues
    eigVals = map(eigVals);
    
    % Get the eigenvalue condition numbers
    condVals = data.condVals;
    
    % Check that the number of values in the real and imaginary files
    % are the same
    if numel(eigVals) ~= numel(condVals)
        error 'Error, eigenvalues and condition number values do not match';
    end
    
    for i = 1:numel(eigVals)
        
        % Get the ith eigenvalue
        z = eigVals(i);
        
        % Get the real and imaginary parts
        xVal = real(z);
        yVal = imag(z);
        
        % Get the condition number
        cVal = condVals(i);
        
        if isInMargin(xVal, yVal, margin)
            
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
            total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
        else
            outsideCount = outsideCount + 1;
        end
    end
    
    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(eigVals)*100);
    
    % Write mesh to a temporary .mat file
    parSave(tmpFilename_count, 1, count);
    parSave(tmpFilename_total, 1, total);

end


% -------------------------------------------------------------------------
function outsideCount = process_cond_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total, map)

    % Sizes of the points
    pointWidth = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;

    % Make the mesh to store the result
    count = uint32(zeros(resolution.height, resolution.width));
    total = zeros(resolution.height, resolution.width);

    % Number of points outside view
    outsideCount = uint64(0);

    % Load the files
    data = parLoad(dataFilename);
    
    % Get the eigenvalues
    eigVals = data.eigVals;
    
    % Map the eigenvalues
    eigVals = map(eigVals);
    
    % Get the eigenvalue condition numbers
    condVals = data.condVals;

    % Check that the number of values in the real and imaginary files
    % are the same
    if numel(eigVals) ~= numel(condVals)
        error 'Error, eigenvalues and condition number values do not match';
    end

    for i = 1:numel(eigVals)

        % Get the ith eigenvalue
        z = eigVals(i);
        
        % Get the real and imaginary parts
        xVal = real(z);
        yVal = imag(z);
        
        % Get the condition number
        cVal = condVals(i);
        
        % 1 =====
        if isInMargin(xVal, yVal, margin)
            
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
            total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
        else
            outsideCount = outsideCount + 1;
        end
        
        % 2 =====
        xVal = xVal * -1;
        if isInMargin(xVal, yVal, margin)
            
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            
            count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
            total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
        else
            outsideCount = outsideCount + 1;
        end
        
        % 3 =====
        yVal = yVal * -1;
        if isInMargin(xVal, yVal, margin)
            
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));
            
            count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
            total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
        else
            outsideCount = outsideCount + 1;
        end
        
        % 4 =====
        xVal = xVal * -1;
        
        if isInMargin(xVal, yVal, margin)
            
            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            
            count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
            total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
        else
            outsideCount = outsideCount + 1;
        end
        
    end
    
    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(eigVals)*100);

    % Write mesh to a temporary .mat file
    parSave(tmpFilename_count, 1, count);
    parSave(tmpFilename_total, 1, total);
    
end


% -------------------------------------------------------------------------
function b = isInMargin(x, y, margin)
    
    b = x > margin.left && x < margin.right && y > margin.bottom && y < margin.top;
    
end


% -------------------------------------------------------------------------
% Use in place of load in a parfor loop
function l = parLoad(fName)
    l = load(fName);
end


% -------------------------------------------------------------------------
% Use in place of save in a parfor loop
function parSave(fname,numvars,varargin)
    for i = 1:numvars
       eval([inputname(i+2),' = varargin{i};']);  
    end
    save('-mat',fname,inputname(3));
    for i = 2:numvars    
        save('-mat',fname,inputname(i+2),'-append');
    end
end