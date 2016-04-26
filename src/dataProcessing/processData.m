% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Apr. 26/2016                                                %
%                                                                         %
% This function will read all .mat files created by the                   %
% generateRandomSample function and sort the eigenvalues onto the complex %
% plane. If the colorByCond input argument is set to Density, then a grid %
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
%   workingDirIn ... (str) The directory where files should be written    %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
%   height ............ Default = 1001 (pixels)                           %
%                       The height (in pixels) of the grid to be used,    %
%                       the width is determined from the margin such that %
%                       each grid point is square                         %
%   margin ............ Default = large enough to fit all the points      %
%                       A struct with keys:                               %
%                            bottom, top, left, right                     %
%                       that indicated the margins for the image          %
%                       if more than one files is used it will only fit   %
%                       all the data in the first file                    %
%   dataFilePrefix .... Default = BHIME                                   %
%                       The prefix for the .mat files that contain the    %
%                       eigenvalues and their condition numbers           %
%   outputFileType .... Default = mat                                     %
%                       Type of file to be output (either mat or txt)     %
%   symmetry .......... Default = false                                   %
%                       If set to true, symmetry across the imaginary     %
%                       axis is used to effectively double the number of  %
%                       points. Symmetry is not used across the real axis %
%                       because for real matrices, the eigenvalues will   %
%                       always appear in conjugate pairs so forcing       %
%                       symmetry would not give any new information.      %
%   numProcessFiles ... Default = All files                               %
%                       The number of files to process                    %
%   map ............... Default = z -> z (no mapping)                     %
%                       Map the eigenvalues by a given function handel    %
%                       The function MUST be vectorized                   %
%   ignoreReal ........ Default = false                                   %
%                       Ignore any points that have zero imaginary part   %
%   ignoreRealTol ..... Default = 1e-8                                    %
%                       The tolerance for how close the imaginary part of %
%                       and eigenvalue is before it is considered to be   %
%                       a real number                                     %
%   colorByCond ....... Default = false                                   %
%                       When set to true, coloring represents the average %
%                       eigenvalue condition number over all the data at  %
%                       each pixel rather than coloring by density        %
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
function [outputFilename, stats] = processData(workingDirIn, options)
    
    tic
    
    narginchk(1, 2)
    
    if nargin < 2
        options = struct();
    end
    
    % Make ProcessData directory if it doesn't exist
    if workingDirIn(end) ~= filesep
        workingDir = [workingDir, filesep];
    else
        workingDir = workingDirIn;
    end
    processDataDir = [workingDir, 'ProcessedData', filesep];
    dataDir = [workingDir, 'Data', filesep];
    mkdir_if_not_exist(processDataDir);
    
    % Process the options
    opts = processOptions(options);
    
    % Get numProcessFiles
    if ~opts.numProcessFilesIsSet
        opts.numProcessFiles = getnumProcessFiles(workingDir, opts.filenamePrefix);
    end
    
    % Get default margin
    if ~opts.marginIsSet
        opts.margin = getDefaultMargin(workingDir, opts.filenamePrefix);
    end
    
    % Get the resolution
    opts.resolution = getResolution(opts.margin, opts.height);
    
    % Output filename
    outputFilename = makeOutputFilename(processDataDir, opts);
    
    % Write readme file
    writeReadMe(processDataDir, outputFilename, opts);
    
    % -------------------------
    
    % Call correct function if colorByCond is true
    if opts.colorByCond
        error('Condition number color is not currently supported.');
        % [mesh, stats] = process_cond();
    else
        [mesh, stats] = process_density(dataDir, opts);
    end
    
    % -------------------------
    
    % Write mesh to a text file
    if strcmp(opts.outputFileType, 'txt')
        dlmwrite([processDataDir, outputFilename], mesh, 'delimiter',' ','precision',15);
    elseif strcmp(opts.outputFileType, 'mat')
        save([processDataDir, outputFilename], 'mesh', 'stats');
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
    totalMesh = uint32(zeros(resolution.height, resolution.width));
    
    % Vector to store number of unique points
    numUniquePts = zeros(1, numProcessFiles);
    
    % Sizes of the points
    pointWidth  = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;
    
    % Process all the data files
    parfor i = 1:numProcessFiles
        
        localMesh = uint32(zeros(resolution.height, resolution.width));
        
        fprintf('File %d of %d\n', i, numProcessFiles);
        
        dataFilename = [dataDir, filenamePrefix, ...
                        '_', num2str(i), '.mat'];
        
        % Load the eigenvalues
        z = parLoad(dataFilename);
        z = z.eigVals(:);
        
        % Map the eigenvalues
        z = opts.map(z);
        
        % tolerance for ignoreReal option
        tol = opts.ignoreRealTol;
        
        % If symmetry, add reflection of values across imaginary axis
        if opts.symmetry
            z = [z, -conj(z)];
        end
        
        % Remove points not in margin
        z = z(isInMargin(real(z), imag(z), margin));
        
        % Ignore real points if ignoreReal option is true
        if opts.ignoreReal
            z = z(abs(imag(z)) > tol);
        end
        
        xVal = uint32(ceil((real(z) - margin.left)/pointWidth));
        yVal = uint32(ceil((imag(z) - margin.bottom)/pointHeight));
        
        idx = uint32(((xVal - 1)*resolution.height + yVal));
        
        % Count number of occurrences of each unique value
        y = sort(idx);
        p = find([true; diff(y)~=0; true]);
        values = y(p(1:end-1));
        instances = diff(p);
        
        % Increment appropriate values
        localMesh(values) = uint32(instances);
        
        totalMesh = totalMesh + localMesh;
        
        %numUniquePts(i) = length(totalMesh(totalMesh ~= 0));
        
    end
    
    stats = struct();
    
end


% ------------------------------------------------------------------- %
% writeReadMe                                                         %
%                                                                     %
% Write a readme file in the processDataDir with information about    %
% when the data was created, what was used to create the data, etc.   %
% ------------------------------------------------------------------- %
function writeReadMe(processDataDir, outputFilename, opts)
    
    % Get the options
    margin          = opts.margin;
    map             = opts.map;
    resolution      = opts.resolution;
    symmetry        = opts.symmetry;
    colorByCond     = opts.colorByCond;
    numProcessFiles = opts.numProcessFiles;
    ignoreReal      = opts.ignoreReal;
    ignoreRealTol   = opts.ignoreRealTol;
    
    % Write a readme file
    file = fopen([processDataDir, 'README.txt'],'a');
    
    % File name
    fprintf(file, [outputFilename, '\n']);
    
    % Date
    fprintf(file, '    Created ............. %s\n', ...
                   datestr(now,'mmmm dd/yyyy HH:MM:SS AM'));
    
    % Margin
    fprintf(file, '    margin .............. bottom: %.5f\n', ...
                   margin.bottom);
    fprintf(file, '                             top: %.5f\n', ...
                   margin.top);
    fprintf(file, '                            left: %.5f\n', ...
                   margin.left);
    fprintf(file, '                           right: %.5f\n', ...
                   margin.right);
                   
    % Map
    fprintf(file, '    map ................. %s\n', func2str(map));
    
    % Resolution
    fprintf(file, '    resolution .......... %d x %d\n', ...
                   resolution.width, resolution.height);
    
    % Symmetry
    fprintf(file, '    symmetry ............ ');
    if symmetry
        fprintf(file, 'true\n');
    else
        fprintf(file, 'false\n');
    end
    
    % colorByCond
    fprintf(file,  '    colorByCond ......... ');
    if colorByCond
        fprintf(file, 'true\n');
    else
        fprintf(file, 'false\n');
    end
    
    % numProcessFiles
    fprintf(file, '    numProcessFiles ..... %d\n', numProcessFiles);
    
    % ignoreReal
    fprintf(file, '    ignoreReal .......... ');
    if ignoreReal
        fprintf(file, 'true\n');
        fprintf(file, '    ignoreRealTol ....... %.5E\n', ...
                       ignoreRealTol);
    else
        fprintf(file, 'false\n');
    end
    
    fprintf(file, '\n\n\n');
    fclose(file);
    
end


% -------------------------------------------------------------------------
% process_density_no_symmetry_tmp
% -------------------------------------------------------------------------
function process_density_tmp(dataFilename, tmpFilename, opts)
    
    resolution = opts.resolution;
    margin     = opts.margin;
    
    % Sizes of the points
    pointWidth  = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;
    
    % Make the mesh to store the result (local to loop)
    mesh = uint32(zeros(resolution.height, resolution.width));
    
    % Load the eigenvalues
    z = parLoad(dataFilename);
    z = z.eigVals(:);
    
    % Map the eigenvalues
    z = opts.map(z);
    
    % tolerance for ignoreReal option
    tol = opts.ignoreRealTol;
    
    % If symmetry, add reflection of values across imaginary axis
    if opts.symmetry
        z = [z, -conj(z)];
    end
    
    z = z(isInMargin(real(z), imag(z), margin));
    
    if opts.ignoreReal
        z = z(abs(imag(z)) > tol);
    end
    
    xVal = uint32(ceil((real(z) - margin.left)/pointWidth));
    yVal = uint32(ceil((imag(z) - margin.bottom)/pointHeight));
    
    idx = uint32(((xVal - 1)*resolution.height + yVal));
    
    % Count number of occurrences of each unique value
    y = sort(idx);
    p = find([true; diff(y)~=0; true]);
    values = y(p(1:end-1));
    instances = diff(p);
    
    % Increment appropriate values
    %mesh(values) = mesh(values) + uint32(instances);
    mesh(values) = uint32(instances);
    
    % Write mesh to a temporary .mat file
    parSave(tmpFilename, 1, mesh);

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
function n = getnumProcessFiles(workingDir, filenamePrefix)
    
    % Data directory
    if workingDir(end) == filesep
        dataDir = [workingDir, 'Data', filesep];
    else
        dataDir = [workingDir, filesep, 'Data', filesep];
    end
    
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
function m = getDefaultMargin(workingDir, dataFilePrefix)
    
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


% ----------------------------------------------------------------------- %
% outputFilename                                                          %
%                                                                         %
% Determine the name of the output file.                                  %
%                                                                         %
% OUTPUT                                                                  %
%   A string of the form                                                  %
%   {dataFilePrefix}-{width}x{height}-{colorBy}-{symmetry}                %
% ----------------------------------------------------------------------- %
function outputFilename = makeOutputFilename(processDataDir, opts)
    
    resolution      = opts.resolution;
    filenamePrefix  = opts.filenamePrefix;
    colorByCond     = opts.colorByCond;
    numProcessFiles = opts.numProcessFiles;
    outputFileType  = opts.outputFileType;
    symmetry        = opts.symmetry;
    
    outPrefix = [filenamePrefix, '-', ...
                 num2str(resolution.width), 'x', ...
                 num2str(resolution.height)];
    if colorByCond
        outPrefix = [outPrefix, '-Cond'];
    end
    
    if symmetry
        outPrefix = [outPrefix, '-sym-', num2str(numProcessFiles)];
    else
        outPrefix = [outPrefix, '-', num2str(numProcessFiles)];
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
















    
    
    
    % ---------------------------------------------------------------------
    %{
    function process_data_files_density()
        
        numProcessFiles = opts.numProcessFiles;
        filenamePrefix  = opts.filenamePrefix;
        
        % Process all the data files
        parfor i = 1:numProcessFiles
            
            fprintf('File %d of %d\n', i, numProcessFiles);
            
            tmpFilename  = [processDataDir, ...
                            'tmp_', num2str(i), '.mat'];
            
            dataFilename = [dataDir, filenamePrefix, ...
                            '_', num2str(i), '.mat'];
            
            % Call function to save the processed data to a temporary file
            process_density_tmp(dataFilename, tmpFilename, opts);
        end
        
    end
    %}
    
    
    % ---------------------------------------------------------------------
    %{
    function [meshs, stats] = total_tmp_files_density()
        
        resolution      = opts.resolution;
        numProcessFiles = opts.numProcessFiles;
        
        % Make the mesh to store the result
        meshs = uint32(zeros(resolution.height, resolution.width));
        
        % Vector to store number of unique points
        numUniquePts = zeros(1, numProcessFiles);
        
        disp('THIS PART');
        
        % Read all temprary files and combine
        parfor i = 1:numProcessFiles
            
            % Get the mesh from the temporary
            f = parLoad([processDataDir, 'tmp_', num2str(i), '.mat']);
            
            % Update cumulative mesh
            meshs = meshs + f.mesh;
            
            % Count number of unique points
            %numUniquePts(i) = length(meshs(meshs ~= 0));
            
            % Delete the temporary file
            %delete([processDataDir, 'tmp_', num2str(i), '.mat']);
            
        end
        
        % -------------------------
        % Fill in the statistics struct
        stats = struct();
        
        % Number of unique points
        stats.numUniquePts = numUniquePts;
        
    end
    %}
    
    % ---------------------------------------------------------------------
    %{
    function [mesh, stats] = process_cond()
        
        % Process the data files into temporary files
        process_data_files_cond();
        
        % Add up all the temporary files
        [mesh, stats] = total_tmp_files_cond();
        
    end
    %}
    %{
    % ---------------------------------------------------------------------
    function process_data_files_cond()
        
        % Process all the data files
        parfor i = 1:numProcessFiles
            
            fprintf('File %d of %d\n', i, numProcessFiles);
            
            tmpFilename_count  = [processDataDir, ...
                                  'tmp_count_', ...
                                  num2str(i), '.mat'];
            tmpFilename_total  = [processDataDir, ...
                                  'tmp_total_', ...
                                  num2str(i), '.mat'];
            
            dataFilename = [dataDir, dataFilePrefix, ...
                            '_', num2str(i), '.mat'];
            
            % Call function to save the processed data to a temporary file
            process_cond_tmp(dataFilename, ...
                             tmpFilename_count, ...
                             tmpFilename_total, ...
                             opts);
        end
        
    end
    
    %}
    %{
    % ---------------------------------------------------------------------
    function [mesh, stats] = total_tmp_files_cond()
        
        % Make the mesh to store the result
        count = uint32(zeros(resolution.height, resolution.width));
        total =        zeros(resolution.height, resolution.width);
        mesh  =        zeros(resolution.height, resolution.width);
        
        % Vector to store number of unique points
        numUniquePts = zeros(1, numProcessFiles);
        
        % Read all temprary files and combine
        for i = 1:numProcessFiles
            
            c = load([processDataDir, ...
                      'tmp_count_', ...
                      num2str(i), '.mat'], 'count');
            t = load([processDataDir, ...
                      'tmp_total_', ...
                      num2str(i), '.mat'], 'total');
            
            count = count + c.count;
            total = total + t.total;
            
            numUniquePts(i) = length(t.total(t.total ~= 0));
            
            delete([processDataDir, 'tmp_count_', num2str(i), '.mat']);
            delete([processDataDir, 'tmp_total_', num2str(i), '.mat']);
            
        end
        
        % Compute the average
        for i=1:numel(mesh)
            if count(i) ~= 0
                mesh(i) = total(i)/double(count(i));
                
                if mesh(i) == 0 && total(i) ~= 0
                    mesh(i) = realmin;
                end
            end
        end
        
        % -------------------------
        % Fill in the statistics struct
        stats = struct();
        
        % Number of unique points
        stats.numUniquePts = numUniquePts;
        
    end
    %}
    
% -------------------------------------------------------------------------
function process_cond_tmp(dataFilename, tmpFilename_count,  tmpFilename_total, opts)
end
%{
    resolution = opts.resolution;
    margin = opts.margin;

    % Sizes of the points
    pointWidth  = (margin.right - margin.left)/resolution.width;
    pointHeight = (margin.top - margin.bottom)/resolution.height;

    % Make the mesh to store the result
    count = uint32(zeros(resolution.height, resolution.width));
    total = zeros(resolution.height, resolution.width);

    % Load the eigenvalues
    data = parLoad(dataFilename);
    z = data.eigVals;

    % Map the eigenvalues
    z = opts.map(z);

    % Get the eigenvalue condition numbers
    condVals = data.condVals;

    % tolerance for ignoreReal option
    tol = opts.ignoreRealTol;

    % If symmetry, add reflection of values across imaginary axis
    if opts.symmetry
        z = [z, -conj(z)];
        condVals = [condVals, condVals];
    end

    valid = isInMargin(real(z), imag(z), margin);

    z = z(valid);
    condVals = condVals(valid);

    if opts.ignoreReal
        valid = abs(imag(z)) > tol;
        z = z(valid);
        condVals = condVals(valid);
    end

    idx = uint32((ceil((real(z) - margin.left)/pointWidth) - 1)*resolution.height + ceil((imag(z) - margin.bottom)/pointHeight));

    % Count number of occurrences of each unique value
    [y, I] = sort(idx(:));
    condVals = condVals(I);
    p = find([true;diff(y)~=0;true]);
    values = y(p(1:end-1));
    instances = diff(p);

    count(values) = count(values) + uint32(instances);
    total(values) = total(values) + condVals(instances)


    if isInMargin(xVal, yVal, margin) && (abs(yVal) > tol || ~ignoreReal)

            xIdx = uint32(ceil((xVal - margin.left)/pointWidth));
            yIdx = uint32(ceil((yVal - margin.bottom)/pointHeight));

            count(yIdx, xIdx) = count(yIdx, xIdx) + 1;
            total(yIdx, xIdx) = total(yIdx, xIdx) + cVal;
        else
            outsideCount = outsideCount + 1;
        end
    end

    % Write mesh to a temporary .mat file
    parSave(tmpFilename_count, 1, count);
    parSave(tmpFilename_total, 1, total);

end
%}