% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Mar. 6/2016                                                 %
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
%   height ....... (int) height of grid (in pixles)                       %
%   margin ....... (struct) {bottom, top, left, right}                    %
%   workingDir ... (str) The directory where files should be written      %
%   colorBy ...... Either 'density' or 'cond'                             %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
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
%                                                                         %
% OUTPUT                                                                  %
%   A string of the name of file that data is written to                  %
%                                                                         %
% TO DO                                                                   %
%   - Separate symmetry into symmetry across real axis and symmetry       %
%     across imaginary axis                                               %
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
function outputFilename = processData(height, margin, workingDir, colorBy, options)
    
    tic
    
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
    
    % Make ProcessData directory if it doesn't exist
    if workingDir(end) == filesep
        processDataDir = [workingDir, 'ProcessedData', filesep];
        dataDir = [workingDir, 'Data', filesep];
    else
        processDataDir = [workingDir, filesep, 'ProcessedData', filesep];
        dataDir = [workingDir, filesep, 'Data', filesep];
    end
    mkdir_if_not_exist(processDataDir)
    
    % Process the options
    opts = processOptions();
    dataFilePrefix = opts.dataFilePrefix;
    outputFileType = opts.outputFileType;
    symmetry = opts.symmetry;
    numFiles = opts.numFiles;
    
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
        mesh = process_density();
    elseif strcmpi(colorBy, 'cond')
        % Call the method for condition number coloring
        mesh = process_cond();
    else
        error 'Invalid value supplied for colorBy argument'
    end

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
    
    
    % ----------------------------------------------------------------------
    function mesh = process_density()
        if symmetry
            mesh = process_density_symmetry();
        else
            mesh = process_density_no_symmetry();
        end
    end
        
    % ----------------------------------------------------------------------
    
    function mesh = process_cond()
        if symmetry
            mesh = process_cond_symmetry();
        else
            mesh = process_cond_no_symmetry();
        end
    end
    
    % ------------------------------------------------------------------
    function mesh = process_density_no_symmetry()
    
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename  = [processDataDir, 'tmp_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat'];
    
            % Call function to save the processed data to a temporary file
            process_density_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename);
    
            toc
    
        end
    
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

    
    % ----------------------------------------------------------------------
    function mesh = process_density_symmetry()
    
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename  = [processDataDir, 'tmp_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat'];
    
            % Call function to save the processed data to a temporary file
            process_density_symmetry_tmp(resolution, margin, dataFilename, tmpFilename)
    
            toc
    
        end
    
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
    
    
    % ======================================================================
    % ======================================================================
    
    function mesh = process_cond_no_symmetry()
    
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename_count  = [processDataDir, 'tmp_count_',num2str(i),'.mat'];
            tmpFilename_total  = [processDataDir, 'tmp_total_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat']
    
            % Call function to save the processed data to a temporary file
            process_cond_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total);
    
            toc
    
        end
    
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
    
    % ----------------------------------------------------------------------
    
    % ----------------------------------------------------------------------
    function mesh = process_cond_symmetry()
    
        % Process all the data
        parfor i = 1:numFiles
    
            tic
    
            disp(['file ', num2str(i), ' of ', num2str(numFiles)]);
    
            tmpFilename_count  = [processDataDir, 'tmp_count_',num2str(i),'.mat'];
            tmpFilename_total  = [processDataDir, 'tmp_total_',num2str(i),'.mat'];
    
            dataFilename = [dataDir, dataFilePrefix, '_', num2str(i), '.mat'];
    
            % Call function to save the processed data to a temporary file
            process_cond_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total);
    
            toc
    
        end
    
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
    
    % ----------------------------------------------------------------------    
    
    % =====================================================================
    % MORE FUNCTION
    % =====================================================================
    
    
    % ------------------------------------------------------------------- %
    % processOption                                                       %
    %                                                                     %
    % Process the options input options struct. If an option is not in    %
    % the options struct the default value is used.                       %
    %                                                                     %
    % INPUT                                                               %
    %   options ... (struct) contains keys corresponding to the options   %
    %                                                                     %
    % OUTPUT                                                              %
    %   A struct opts with keys                                           %
    %       dataFilePredix                                                %
    %       outputFileType                                                %
    %       symmetry                                                      %
    %       numFiles                                                      %
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
        
        optionNames = struct('dataFilePrefix', 'dataFilePrefix', ...
                             'outputFileType', 'outputFileType', ...
                             'symmetry', 'symmetry', ...
                             'numFiles', 'numFiles');
        
        if ~all(ismember(fnames, fieldnames(optionNames)))
            error('processData:InvalidOption',  ...
                  'Invalid option provided');
        end
        
        % Default values
        dataFilePrefix = 'BHIME';
        outputFileType = 'mat';
        symmetry = false;
        numFiles = getNumFiles();
        
        % dataFilenamePredix
        if isfield(options, 'dataFilePrefix')
            dataFilePrefix = options.dataFilePredix;
        end
        
        % outputFileType
        if isfield(options, 'outputFileType')
            outputFileType = options.outputFileType;
        end
        
        % symmetry
        if isfield(options, 'symmetry')
            symmetry = options.symmetry;
        end
        
        % numFiles
        if isfield(options, 'numFiles')
            numFiles = options.numFiles;
        end
        
        opts = struct('dataFilePrefix', dataFilePrefix, ...
                      'outputFileType', outputFileType, ...
                      'symmetry', symmetry, ...
                      'numFiles', numFiles);
        
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
            error 'Bottom margin must be less than top margin'
        end
        if margin.left >= margin.right
            error 'Left margin must be less than top margin'
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






% ------------------------------------------------------------------
function process_density_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename)

    % Sizes of the points
    pointWidth = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;

    % Make the mesh to store the result (local to loop)
    mesh = uint32(zeros(resolution.height, resolution.width));

    % Number of points outside view
    outsideCount = uint64(0);

    % Load the files
    data = parLoad(dataFilename);

    for i = 1:numel(data.eig)

        xVal = real(data.eig(i));
        yVal = imag(data.eig(i));
        
        if xVal > margin.left && xVal <= margin.right && yVal > margin.bottom && yVal <= margin.top

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx <= resolution.width && yIdx > 0 && yIdx <= resolution.height
                mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
            end
        else
            outsideCount = outsideCount + 1;
        end
    end

    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(data.eig)*100);

    % Write mesh to a temporary .mat file
    parSave(tmpFilename, 1, mesh);

end


% ----------------------------------------------------------------------
function process_density_symmetry_tmp(resolution, margin, dataFilename, tmpFilename)

    % Sizes of the points
    pointWidth = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;

    % Make the mesh to store the result (local to loop)
    mesh = uint32(zeros(resolution.height, resolution.width));

    % Number of points outside view
    outsideCount = uint64(0);

    % Load the files
    data = parLoad(dataFilename);

    for i = 1:numel(data.eig)

        xVal = real(data.eig(i));
        yVal = imag(data.eig(i));

        
        if xVal > margin.left && xVal <= margin.right && yVal > margin.bottom && yVal <= margin.top

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx <= resolution.width && yIdx > 0 && yIdx <= resolution.height
                mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
            end

            % 2 =====
            xVal = xVal * -1;

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx < resolution.width && yIdx > 0 && yIdx < resolution.height
                mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
            end

            % 3 =====
            yVal = yVal * -1;

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx < resolution.width && yIdx > 0 && yIdx < resolution.height
                mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
            end

            % 4 =====
            xVal = xVal * -1;

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx < resolution.width && yIdx > 0 && yIdx < resolution.height
                mesh(yIdx, xIdx) = mesh(yIdx, xIdx) + 1;
            end

        else
            outsideCount = outsideCount + 1;
        end
    end

    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(data.eig)*100);

    % Write mesh to a temporary .mat file
    parSave(tmpFilename, 1, mesh);

end


function process_cond_no_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total)

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

    % Check that the number of values in the real and imaginary files
    % are the same
    if numel(data.eig) ~= numel(data.cond)
        error 'Error, eigenvalues and condition number values do not match';
    end

    for i = 1:numel(data.eig)

        xVal = real(data.eig(i));
        yVal = imag(data.eig(i));
        cVal = data.cond(i);

        if xVal > margin.left && xVal <= margin.right && yVal > margin.bottom && yVal <= margin.top

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx <= resolution.width && yIdx > 0 && yIdx <= resolution.height
                count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
                total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
            end
        else
            outsideCount = outsideCount + 1;
        end
    end

    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(data.eig)*100);
    
    % Write mesh to a temporary .mat file
    parSave(tmpFilename_count, 1, count);
    parSave(tmpFilename_total, 1, total);

end

function process_cond_symmetry_tmp(resolution, margin, dataFilename, tmpFilename_count, tmpFilename_total)

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

    % Check that the number of values in the real and imaginary files
    % are the same
    if numel(data.eig) ~= numel(data.cond)
        error 'Error, eigenvalues and condition number values do not match';
    end

    for i = 1:numel(data.eig)

        xVal = real(data.eig(i));
        yVal = imag(data.eig(i));
        cVal = data.cond(i);


        if xVal > margin.left && xVal <= margin.right && yVal > margin.bottom && yVal <= margin.top

            % 1 =====
            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx <= resolution.width && yIdx > 0 && yIdx <= resolution.height
                count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
                total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
            end

            % 2 =====
            xVal = xVal * -1;

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx < resolution.width && yIdx > 0 && yIdx < resolution.height
                count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
                total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
            end

            % 3 =====
            yVal = yVal * -1;

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx < resolution.width && yIdx > 0 && yIdx < resolution.height
                count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
                total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
            end

            % 4 =====
            xVal = xVal * -1;

            xIdx = uint32(floor((xVal - margin.left)/pointWidth));
            yIdx = uint32(floor((yVal - margin.bottom)/pointHeight));
            yIdx = resolution.height - yIdx + 1;

            if xIdx > 0 && xIdx < resolution.width && yIdx > 0 && yIdx < resolution.height
                count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
                total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
            end

        else
            outsideCount = outsideCount + 1;
        end
    end

    % Print number of points outside view
    fprintf('%.4f%% of points outside margins\n', outsideCount/numel(data.eig)*100);

    % Write mesh to a temporary .mat file
    parSave(tmpFilename_count, 1, count);
    parSave(tmpFilename_total, 1, total);
    
end



% Use in place of load in a parfor loop
function l = parLoad(fName)
    l = load(fName);
end


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