% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Feb. 19/2018                                                %
%                                                                         %
% This function will generate .mat files containing the for a sample of   %
% matrices. The matrices are sampled from a user provided function.       %
%   - A new directory 'Data' will be created in the input workingDir, or  %
%     if the option overrideDataDir is set, this will be created and used %
%     to store the files created by this function.                        %
%   - A readme.txt file will be automatically generated in the data       %
%     directory.                                                          %
%                                                                         %
% INPUT                                                                   %
%   generator ............... A function handle that takes no input and   %
%                             returns  a random square matrix, must       %
%                             always return a matrix of the same size.    %
%   workingDirIn ............ The directory where files should be written %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
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
%                             Control how many matrices eigenvalues are   %
%                             in each data file.                          %
%   numDataFiles ............ Default = 1                                 %
%                             Set this option to a positive integer if    %
%                             you would like to generate multiple files   %
%                             with data where each file contains the      %
%                             eigenvalues for matricesPerFile random      %
%                             matrices.                                   %
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
    
    if ~opts.startFileIndexIsSet
        startFileIndex = getStartFileIndex(dataDir, filenamePrefix);
    end
    
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
        
        eigVals = computeEig(generator, matrixSize, opts);
        
        % Save the data
        generatorStr = func2str(generator);
        filename = [dataDir, filenamePrefix, '_', num2str(k), '.mat'];
        
        save(filename, 'eigVals', ...
                       'matrixSize', ...
                       'generatorStr', ...
                       '-v7.3');
        
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
    innerLoopSize = max(min(1e4, floor(opts.matricesPerFile/(10*feature('numCores')))),1);
    
    % Outer loop size
    outerLoopSize = opts.matricesPerFile/innerLoopSize;
    
    % Vector to store the computed eigenvalues
    eigVals = zeros(innerLoopSize*matrixSize, outerLoopSize);
    
    parfor i=1:outerLoopSize
        
        eigValsLocal = zeros(innerLoopSize, matrixSize);
        
        for j=1:innerLoopSize
            try
                eigValsLocal(j, :) = eig(generator());
            end
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