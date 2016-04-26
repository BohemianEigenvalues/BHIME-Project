% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Apr. 26/2016                                                %
%                                                                         %
% This function will generate .mat files containing the eigenvalues and   %
% their respective condition numbers for a sample of matrices. The        %
% matrices are sampled from a user provided function.                     %
%   - A readme file will be automatically generated when this function    %
%     is called                                                           %
%   - A new directory 'Data' in the input workingDir will be created      %
%   - Files will be named BHIME_i.mat where i starts at 1 (you can        %
%     customize the prefix by using the filenamePrefix option)            %
%                                                                         %
% INPUT                                                                   %
%   generator ......... A function handle that takes no input and returns %
%                       a random square matrix, must always return a      %
%                       matrix of the same size                           %
%   workingDir ........ The directory where files should be written       %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
%   filenamePrefix .... Default = BHIME                                   %
%                       The prefix for the .mat files that contain the    %
%                       eigenvalues and their condition numbers           %
%   startFileIndex .... Default = one more than the last file index if    %
%                                 there are any previously generated      %
%                                 files, 1 if there are no previously     %
%                                 generated files.                        %
%                       If you want your files to be indexed with a       %
%                       custom starting value set this value to a         %
%                       positive integer value.                           %
%   numDataFiles ...... Default = 1                                       %
%                       Set this value if you would like to generate      %
%                       multiple data files                               %
%   matricesPerFile ... Default = floor(1e6/matrixSize)                   %
%                       Set this value to control how many matrices       %
%                       eigenvalues will be in each file, for larger      %
%                       matrices a smaller value will be needed. You can  %
%                       determine a good value based on the size of       %
%                       matrix you're using and the amount of memory your %
%                       computer has.                                     %
%                          file size = memory usage =                     %
%                             64*matrixSize*matricesPerFile bits          %
%                             (Eigenvalues and condition numbers are      %
%                             stored in single precision)                 %
%   computeCond ....... Default = true                                    %
%                       Set this value to false to only compute the       %
%                       eigenvalues. Should result in a speed improvement %
%                       of roughly 6 times faster.                        %
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
    
    % Make Data directory if it doesn't exist
    if workingDirIn(end) ~= filesep
        workingDir = [workingDirIn, filesep];
    else
        workingDir = workingDirIn;
    end
    dataDir = [workingDir, 'Data', filesep];
    mkdir_if_not_exist(dataDir);
    
    % Determine matrix size
    matrixSize = max(size(generator()));
    
    % Process the options
    opts = processOptions(options);
    
    startFileIndex  = opts.startFileIndex;
    numFiles        = opts.numDataFiles;
    matricesPerFile = opts.matricesPerFile;
    filenamePrefix  = opts.filenamePrefix;
    computeCond     = opts.computeCond;
    
    if ~opts.startFileIndexIsSet
        startFileIndex = getStartFileIndex(dataDir, filenamePrefix)
    end
    if ~opts.matricesPerFileIsSet
        matricesPerFile = floor(1e6/matrixSize);
    end
    
    % Check for a parameters.mat file
    checkForParametersFile(dataDir, generator, matrixSize, opts);
    
    % Write a readme file
    writeReadMe(dataDir, generator, matrixSize, opts)
    
    % Set random seed based on the current time 
    s = RandStream('mt19937ar', 'Seed', mod(int64(now*1e6), 2^32));
    RandStream.setGlobalStream(s);
    
    % Master timer
    tic;
    
    % Last index for file to be written
    endFileIndex = startFileIndex + numFiles - 1;
    
    for k=startFileIndex:endFileIndex
        
        % Timer 2
        t2 = tic;
        
        fprintf('\n%d of %d\n', k, endFileIndex);
        
        if computeCond
            [eigVals, condVals] = computeEigAndCond();
        else
            eigVals = computeEig();
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
    averageTime = toc/numFiles;
    fprintf('\nAverage loop time is %.3f seconds\n', averageTime);
    
    
    % =====================================================================
    % FUNCTIONS
    % =====================================================================
    
    
    % ------------------------------------------------------------------- %
    % computeEigAndCond                                                   %
    %                                                                     %
    % Compute both the eigenvalues and their condition numbers            %
    % ------------------------------------------------------------------- %
    function [eigVals, condVals] = computeEigAndCond()
        
        % Vectors to store eigenvalues and condition numbers
        eigVals  = single(zeros(matricesPerFile, matrixSize));
        condVals = single(zeros(matricesPerFile, matrixSize));
        
        parfor i=1:matricesPerFile
            
            % Generator a random matrix
            A = generator();
            
            % Compute eigenvalues and condition numbers w.r.t. eigenvalues
            [~, D, s] = condeig(A);
            
            eigVals(i,:)  = single(diag(D));
            condVals(i,:) = single(s);
            
        end
        
    end
    
    
    % ------------------------------------------------------------------- %
    % computeEig                                                          %
    %                                                                     %
    % Compute only the eigenvalues                                        %
    % ------------------------------------------------------------------- %
    function eigVals = computeEig()
        
        % Vectors to store eigenvalues and condition numbers
        eigVals = single(zeros(matricesPerFile, matrixSize));
        
        parfor i=1:matricesPerFile
            
            eigVals(i, :) = single(eig(generator()));
            
        end
        
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
% values of several parameters specific to the data:                      %
%   generator ......... Function handle used to generate random matrices  %
%   matrixSize ........ Size of the matrices that are generated           %
%   filenamePrefix .... Prefix for the data files                         %
%   matricesPerFile ... Number of matrices eigenvalues in each file       %
%   computeCond ....... Were condition numbers computed?                  %
%                                                                         %
% INPUT                                                                   %
%   dataDir ...... Directory with data                                    %
%   generator .... Matrix generator function handle                       %
%   matrixSize ... Size of the matrices that are generated                %
%   opts ......... The option struct                                      %
% ----------------------------------------------------------------------- %
function writeParametersFile(dataDir, generator, matrixSize, opts)
    
    filenamePrefix = opts.filenamePrefix;
    matricesPerFile = opts.matricesPerFile;
    computeCond = opts.computeCond;
    
    save([dataDir, 'parameters.mat'], ...
         'generator', ...
         'matrixSize', ...
         'filenamePrefix', ...
         'matricesPerFile', ...
         'computeCond');
    
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
    
    % Check filenamePrefix
    if ~strcmp(f.filenamePrefix, opts.filenamePrefix)
        fprintf('filenamePrefix does not match\n');
        b = false;
    end
    
    % Check matricesPerFile
    if f.matricesPerFile ~= opts.matricesPerFile
        fprintf('matricesPerFile does not match\n');
        b = false;
    end
    
    % Check computeCond
    if f.computeCond ~= opts.computeCond
        fprintf('computeCond does not match\n');
        b = false;
    end
    
    % Check generator
    if ~strcmp(f.generator, generator)
        fprintf('generator does not match\n');
        b = false;
    end
    
    % Check matrixSize
    if f.matrixSize ~= matrixSize
        fprintf('matrixSize does not match\n');
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
    fprintf(file, 'Last Updated ...... %s\n', ...
            datestr(now,'mmmm dd/yyyy HH:MM:SS AM'));
    
    % matricesPerFile
    fprintf(file, 'matricesPerFile ... %d\n', opts.matricesPerFile);
    
    % computeCond
    fprintf(file, 'computeCond ....... ');
    if opts.computeCond
        fprintf(file, 'true\n');
    else
        fprintf(file, 'false\n');
    end
    
    % generator
    fprintf(file, ['generator ......... ', func2str(generator), '\n']);
    
    % matrixSize
    fprintf(file, 'matrixSize ........ %d\n', matrixSize);
    
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