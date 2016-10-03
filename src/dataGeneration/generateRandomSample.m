% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Oct. 3/2016                                                 %
%                                                                         %
% This function will generate .mat files containing the eigenvalues and   %
% their condition numbers (if option computeCond is true) for a sample of %
% matrices. The matrices are sampled from a user provided function.       %
%   - A new directory 'Data' will be created in the input workingDir, or  %
%     if the option overrideDataDir is set, this will be created and used %
%     to store the files created by this function.                        %
%   - A readme.txt file will be automatically generated in the data       %
%     directory.                                                          %
%   - A parameters.mat file will be automatically generated in the data   %
%     directory. This will contain information about the parameter        %
%     (option) values used when calling the function and ensures          %
%     consistency when the more data is added to an existing directory    %
%                                                                         %
% INPUT                                                                   %
%   generator ............... A function handle that takes no input and   %
%                             returns  a random square matrix, must       %
%                             always return a matrix of the same size.    %
%   workingDirIn ............ The directory where files should be written %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
%   computeCond ............. Default = false                             %
%                             When set to true, this function will        %
%                             compute and store the condition numbers of  %
%                             the eigenvalues. This will take roughly 6   %
%                             times longer to execute when set to true.   %
%   dataPrecision ........... Default = 'single'                          %
%                             This can be set to either 'single' or       %
%                             'double'. This option does not affect the   %
%                             precision the eigenvalues are computed in.  %
%                             Eigenvalues are always computed in double   %
%                             precision. By setting this to 'single' the  %
%                             eigenvalues will be stored in single        %
%                             precision. That is, they will be  computed  %
%                             in double precision and cast to single      %
%                             precision for storage.                      %
%   filenamePrefix .......... Default = 'BHIME'                           %
%                             The name that will be used when naming the  %
%                             data files. The names of the data files     %
%                             take the form: filenamePrefix + '_' + i     %
%                             where i is a positive integer.              %
%   ignoreRealData .......... Default = false                             %
%                             When this option is true, eigenvalues where %
%                             the  imaginary part is within ignoreRealTol %
%                             of zero will not be stored in the data      %
%                             file. By setting this to true you can       %
%                             substantially reduce the size of the data   %
%                             files.                                      %
%   ignoreRealTol ........... Default = 1e-10                             %
%                             A complex number z will be considered a     %
%                             real value, and therefore not stored in the %
%                             output data file (if ignoreRealData is      %
%                             true) if abs(Im(z)) < ignoreRealTol.        %
%   matricesPerFile ......... Default = floor(1e6/matrixSize)             %
%                             Control how many matrices eigenvalues (and  %
%                             condition numbers) are in each data file.   %
%   numDataFiles ............ Default = 1                                 %
%                             Set this option to a positive integer if    %
%                             you would like to generate multiple files   %
%                             with data where each file contains the      %
%                             eigenvalues and their condition numbers for %
%                             matricesPerFile random matrices.            %
%   overrideDataDir ......... Default = [empty string] (not set)          %
%                             Set this option to a non-empty string       %
%                             indicating a directory to write the data    %
%                             files to. If not set they will be written   %
%                             to the directory workingDir/Data/. If the   %
%                             set directory does not exist it will be     %
%                             created.                                    %
%   startFileIndex .......... Default = 1 greater than the  highest index %
%                                       of the files in the data          %
%                                       directory or 1 if no files have   %
%                                       been written.                     %
%                             Only use this if you have already computed  %
%                             data and would like to compute more.        %
%   storeDataWithSymmetry ... Default = false                             %
%                             When this option is set to true, any values %
%                             not in the upper right quadrant of the      %
%                             complex plane will not be stored in the     %
%                             output data files. That is, only values     %
%                             where Im(z) >= 0 and Re(z) >= 0 are stored. %
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
function generateRandomSample(generator, workingDirIn, options)
    
    narginchk(2,3);
    
    if nargin < 3
        options = struct();
    end
    
    % Check that generator is a function handle
    if ~isa(generator, 'function_handle')
        error('generator must be a function handle');
    end
    
    % Process the options
    opts = processOptions(options);
    
    % Make sure the working directory ends in the filesep
    if workingDirIn(end) ~= filesep
        workingDir = [workingDirIn, filesep];
    else
        workingDir = workingDirIn;
    end
    
    % Set and make the data directory
    if opts.overrideDataDirIsSet
        dataDir = opts.overrideDataDir;
    else
        dataDir = [workingDir, 'Data', filesep];
    end
    mkdir_if_not_exist(dataDir);
    
    % Determine matrix size
    matrixSize = max(size(generator()));
    
    startFileIndex  = opts.startFileIndex;
    numDataFiles    = opts.numDataFiles;
    filenamePrefix  = opts.filenamePrefix;
    computeCond     = opts.computeCond;
    
    if ~opts.startFileIndexIsSet
        startFileIndex = getStartFileIndex(dataDir, filenamePrefix);
    end
    
    % Check for a parameters.mat file
    checkForParametersFile(dataDir, generator, matrixSize, opts);
    
    % Write a readme file
    writeReadMe(dataDir, generator, matrixSize, opts);
    
    % Set random seed based on the current time 
    s = RandStream('mt19937ar', 'Seed', mod(int64(now*1e6), 2^32));
    RandStream.setGlobalStream(s);
    
    % Master timer
    tic;
    
    % Last index for file to be written
    endFileIndex = startFileIndex + numDataFiles - 1;
    
    for k=startFileIndex:endFileIndex
        
        % Timer 2
        t2 = tic;
        
        fprintf('\n%d of %d\n', k, endFileIndex);
        
        if computeCond
            [eigVals, condVals] = computeData(generator, matrixSize, opts);
        else
            eigVals = computeData(generator, matrixSize, opts);
        end
        
        % ---------------
        % Save the data
        generatorStr = func2str(generator);
        filename = [dataDir, filenamePrefix, '_', num2str(k), '.mat'];
        
        if computeCond
            save(filename, 'eigVals', ...
                           'condVals', ...
                           'matrixSize', ...
                           'generatorStr');
        else
            save(filename, 'eigVals', ...
                           'matrixSize', ...
                           'generatorStr');
        end
        % ---------------
        
        % Timing
        tElapsed = toc(t2);  % timer 2
        fprintf('Loop time is %.3f seconds\n', tElapsed);
        
    end
    
    % Timing
    averageTime = toc/numDataFiles;
    fprintf('\nAverage loop time is %.3f seconds\n', averageTime);
    
end


% =====================================================================
% FUNCTIONS
% =====================================================================


% ------------------------------------------------------------------- %
% computeData                                                         %
%                                                                     %
% Compute the eigenvalues and optionally their condition numbers      %
% ------------------------------------------------------------------- %
function [eigVals, condVals] = computeData(generator, matrixSize, opts)
    
    if nargout == 1
        eigVals = computeEig(generator, matrixSize, opts);
    else
        [eigVals, condVals] = computeEigAndCond(generator, matrixSize, opts);
    end
    
end


% ------------------------------------------------------------------- %
% computeEigAndCond                                                   %
%                                                                     %
% Compute both the eigenvalues and their condition numbers            %
% ------------------------------------------------------------------- %
function [eigVals, condVals] = computeEigAndCond(generator, matrixSize, opts)
    
    % Size of the inner loop, need to find optimal value for this
    % parameter based on the matricesPerFile variable
    % Maybe matricesPerFile/100?
    % Maybe matricesPerFile/(numCores*10)?
    innerLoopSize = 1e4;
    
    % Outer loop size
    outerLoopSize = opts.matricesPerFile/innerLoopSize;
    
    % Vector to store the computed eigenvalues
    eigVals  = zeros(innerLoopSize*matrixSize, outerLoopSize);
    condVals = zeros(innerLoopSize*matrixSize, outerLoopSize);
    
    parfor i=1:outerLoopSize
        
        eigValsLocal  = zeros(innerLoopSize, matrixSize);
        condValsLocal = zeros(innerLoopSize, matrixSize);
        
        for j=1:innerLoopSize
            
            A = generator();
            
            % Compute eigenvalues and condition numbers w.r.t. eigenvalues
            [~, D, s] = condeig(A);
            
            eigValsLocal(i,:)  = diag(D);
            condValsLocal(i,:) = s;
            
        end
        
        eigVals(:, i)  = eigValsLocal(:);
        condVals(:, i) = condValsLocal(:)
        
    end
    
    % Cast to single in dataPrecision option is set to single
    if strcmp(opts.dataPrecision, 'single')
        eigVals  = single(eigVals);
        condVals = single(condVals);
    end
    
    % Convert to a vector
    eigVals  = eigVals(:);
    condVals = condVals(:);
    
    % ignoreRealData option
    if opts.ignoreRealData
        valid = abs(imag(eigVals)) > opts.ignoreRealTol;
        eigVals  = eigVals(valid);
        condVals = condVals(valid);
        
    end
    
    % storeDataWithSymmetry option
    if opts.storeDataWithSymmetry
        valid = imag(eigVals) >= 0 & real(eigVals) >= 0;
        eigVals  = eigVals(valid);
        condVals = condVals(valid);
    end
    
end


% ------------------------------------------------------------------- %
% computeEig                                                          %
%                                                                     %
% Compute only the eigenvalues                                        %
% Looks a little weird to optimize the speed of using the parfor loop %
% ------------------------------------------------------------------- %
function eigVals = computeEig(generator, matrixSize, opts)
    
    % Size of the inner loop, need to find optimal value for this
    % parameter based on the matricesPerFile variable
    % Maybe matricesPerFile/100?
    % Maybe matricesPerFile/(numCores*10)?
    innerLoopSize = 1e4;
    
    % Outer loop size
    outerLoopSize = opts.matricesPerFile/innerLoopSize;
    
    % Vector to store the computed eigenvalues
    eigVals = zeros(innerLoopSize*matrixSize, outerLoopSize);
    
    parfor i=1:outerLoopSize
        
        eigValsLocal = zeros(innerLoopSize, matrixSize);
        
        for j=1:innerLoopSize
            eigValsLocal(j, :) = eig(generator());
        end
        
        eigVals(:, i) = eigValsLocal(:);
        
    end
    
    % Cast to single in dataPrecision option is set to single
    if strcmp(opts.dataPrecision, 'single')
        eigVals = single(eigVals);
    end
    
    % Convert to a vector
    eigVals = eigVals(:);
    
    % ignoreRealData option
    if opts.ignoreRealData
        eigVals = eigVals(abs(imag(eigVals)) > opts.ignoreRealTol);
    end
    
    % storeDataWithSymmetry option
    if opts.storeDataWithSymmetry
        eigVals = eigVals(imag(eigVals) >= 0 & real(eigVals) >= 0);
    end
    
end


% ----------------------------------------------------------------------- %
% checkForParametersFile                                                  %
%                                                                         %
% Check if a parameters file exists, if it doesn't, create on. If it does %
% verify that the parameters for this call of the function match what is  %
% in the parameters.mat file. If they don't match give the user the       %
% option to abort or overwrite all previous data.                         %
%                                                                         %
% INPUT                                                                   %
%   dataDir ...... Directory with data                                    %
%   generator .... Matrix generator function handle                       %
%   matrixSize ... Size of the matrices that are generated                %
%   opts ......... The option struct                                      %
% ----------------------------------------------------------------------- %
function checkForParametersFile(dataDir, generator, matrixSize, opts)
    
    % Check if parameters.mat file exists
    if exist([dataDir, 'parameters.mat']) == 2
        
        % Verify all parameters match
        match = checkParametersFile(dataDir, generator, matrixSize, opts);
        
        if ~match
            str = '';
            while ~(strcmp(str, 'yes') || strcmp(str, 'no'))
                str = input('Would you like overwrite all old data(yes/no): ', 's');
            end
            
            if strcmp(str, 'yes')
                % Delete all old data files
                delete([dataDir, '*.mat']);
                delete([dataDir, '*.txt']); 
               
                % Write a paramters.mat file
                writeParametersFile(dataDir, generator, matrixSize, opts);
            else
                error('Cannot create data with different parameters.');
            end
            
        end
        
    else
        % Write a paramters.mat file
        writeParametersFile(dataDir, generator, matrixSize, opts);
    end
    
end


% ----------------------------------------------------------------------- %
% writeParametersFile                                                     %
%                                                                         %
% Write a file parameters.mat in the data directory that contains the     %
% values of several parameters specific to the data.                      %
%                                                                         %
% INPUT                                                                   %
%   dataDir ...... Directory with data                                    %
%   generator .... Matrix generator function handle                       %
%   matrixSize ... Size of the matrices that are generated                %
%   opts ......... The option struct                                      %
% ----------------------------------------------------------------------- %
function writeParametersFile(dataDir, generator, matrixSize, opts)
    
    computeCond           = opts.computeCond;
    dataPrecision         = opts.dataPrecision;
    filenamePrefix        = opts.filenamePrefix;
    ignoreRealData        = opts.ignoreRealData;
    ignoreRealTol         = opts.ignoreRealTol;
    matricesPerFile       = opts.matricesPerFile;
    storeDataWithSymmetry = opts.storeDataWithSymmetry;
    
    save([dataDir, 'parameters.mat'], ...
         'generator', ...
         'matrixSize', ...
         'computeCond', ...
         'dataPrecision', ...
         'filenamePrefix', ...
         'ignoreRealData', ...
         'ignoreRealTol', ...
         'matricesPerFile', ...
         'storeDataWithSymmetry');
    
end


% ----------------------------------------------------------------------- %
% checkParametersFile                                                     %
%                                                                         %
% Check if the values in the paramters.mat file match the values used in  %
% this call to the generateRandomSample function.                         %
%                                                                         %
% INPUT                                                                   %
%   dataDir ...... Directory with data                                    %
%   generator .... Matrix generator function handle                       %
%   matrixSize ... Size of the matrices that are generated                %
%   opts ......... The option struct                                      %
% ----------------------------------------------------------------------- %
function b = checkParametersFile(dataDir, generator, matrixSize, opts)
    
    % Load the file
    f = load([dataDir, 'parameters.mat']);
    
    b = true;
    
    % Check computeCond
    if f.computeCond ~= opts.computeCond
        fprintf('computeCond option does not match\n');
        b = false;
    end
    
    % Check dataPrecision
    if f.dataPrecision ~= opts.dataPrecision
        fprintf('dataPrecision option does not match\n');
        b = false;
    end
    
    % Check filenamePrefix
    if ~strcmp(f.filenamePrefix, opts.filenamePrefix)
        fprintf('filenamePrefix option does not match\n');
        b = false;
    end
    
    % Check ignoreRealData
    if f.ignoreRealData ~= opts.ignoreRealData
        fprintf('ignoreRealData option does not match\n');
        b = false;
    end
    
    % Check ignoreRealTol
    if f.ignoreRealTol ~= opts.ignoreRealTol
        fprintf('ignoreRealTol option does not match\n');
        b = false;
    end
    
    % Check matricesPerFile
    if f.matricesPerFile ~= opts.matricesPerFile
        fprintf('matricesPerFile option does not match\n');
        b = false;
    end
    
    % Check storeDataWithSymmetry
    if f.storeDataWithSymmetry ~= opts.storeDataWithSymmetry
        fprintf('storeDataWithSymmetry option does not match\n');
        b = false;
    end
    
    % Check matrixSize
    if f.matrixSize ~= matrixSize
        fprintf('matrixSize does not match\n');
        b = false;
    end
    
    % Check generator
    if ~strcmp(func2str(f.generator), func2str(generator))
        fprintf('generator does not match\n');
        b = false;
    end
    
end


% ----------------------------------------------------------------------- %
% writeReadMe                                                             %
%                                                                         %
% Write a readme file in the dataDir with information about when the data %
% was created, what was used to create the data, etc.                     %
%                                                                         %
% INPUT                                                                   %
%   dataDir ...... Directory with data                                    %
%   generator .... Matrix generator function handle                       %
%   matrixSize ... Size of the matrices that are generated                %
%   opts ......... The option struct                                      %
% ----------------------------------------------------------------------- %
function writeReadMe(dataDir, generator, matrixSize, opts)
    
    % Write a readme file
    file = fopen([dataDir, 'README.txt'], 'w');
    
    % Date
    fprintf(file, 'Last Updated ............ %s\n', ...
            datestr(now,'mmmm dd/yyyy HH:MM:SS AM'));
    
    % computeCond
    fprintf(file, 'computeCond ............. ');
    if opts.computeCond
        fprintf(file, 'true\n');
    else
        fprintf(file, 'false\n');
    end
    
    % dataPrecision
    fprintf(file, 'dataPrecision ........... %s\n', opts.dataPrecision);
    
    % ignoreRealData
    fprintf(file, 'ignoreRealData .......... ');
    if opts.ignoreRealData
        fprintf(file, 'true\n');
        % ignoreRealTol
        fprintf(file, 'ignoreRealTol ........... %.5E\n', ...
                      opts.ignoreRealTol);
    else
        fprintf(file, 'false\n');
    end
    
    % matricesPerFile
    fprintf(file, 'matricesPerFile ......... %d\n', opts.matricesPerFile);
    
    % storeDataWithSymmetry
    fprintf(file, 'storeDataWithSymmetry ... ');
    if opts.storeDataWithSymmetry
        fprintf(file, 'true\n');
    else
        fprintf(file, 'false\n');
    end
    
    % generator
    fprintf(file, 'generator ............... %s\n', func2str(generator));
    
    % matrixSize
    fprintf(file, 'matrixSize .............. %d\n', matrixSize);
    
    fclose(file);
end


% ----------------------------------------------------------------------- %
% getStartFileIndex                                                       %
%                                                                         %
% Determine the index of the last file created.                           %
%                                                                         %
% INPUT                                                                   %
%   dataDir .......... Directory containing data files                    %
%   filenamePrefix ... Prefix for data files                              %
%                                                                         %
% OUTPUT                                                                  %
%   Integer, one more than the last file in the data directory with the   %
%   filenamePrefix prefix.                                                %
% ----------------------------------------------------------------------- %
function imax = getStartFileIndex(dataDir, filenamePrefix)
    
    folderInfo = dir(dataDir);
    
    imax = 0;
    
    fnpLen = numel(filenamePrefix);
    
    for k=1:length(folderInfo)
        
        filename = folderInfo(k).name;
        [~, name, ext] = fileparts(filename);
        
        if ~strcmp(ext, '.mat')
            continue
        end
        if numel(name) < fnpLen + 2
            continue
        end
        if ~strcmp(name(1:fnpLen+1), [filenamePrefix, '_'])
            continue
        end
        
        i = str2double(name(fnpLen+2:end));
        
        imax = max(imax, i);
    end
    
    % Increment by 1
    imax = imax + 1;
    
end