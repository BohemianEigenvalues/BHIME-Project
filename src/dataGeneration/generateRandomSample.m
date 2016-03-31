% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Mar. 7/2016                                                 %
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
%   numFiles .......... Default = 1                                       %
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
function generateRandomSample(generator, workingDir, options)

    % Check number of input arguments
    if nargin > 3
        error('generateRandomSample:TooManyInputs', ...
              'requires at most 3 input arguments');
    elseif nargin < 2
        error('generateRandomSample:NotEnoughInputs', ...
              'requires at least 2 input arguments');
    elseif nargin == 2
        options = struct();
    end
    
    % Check that generator is a function handle
    if ~isa(generator, 'function_handle')
        error('generator must be a function handle');
    end
    
    % Make Data directory if it doesn't exist
    if workingDir(end) == filesep
        dataDir = [workingDir, 'Data', filesep];
    else
        dataDir = [workingDir, filesep, 'Data', filesep];
    end
    mkdir_if_not_exist(dataDir);
    
    % Determine matrix size
    matrixSize = max(size(generator()));
    
    % Process the options
    opts = processOptions();
    startFileIndex = opts.startFileIndex;
    numFiles = opts.numFiles;
    matricesPerFile = opts.matricesPerFile;
    filenamePrefix = opts.filenamePrefix;
    
    % Write a readme file
    writeReadMe();
    
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
        
        % Vectors to store eigenvalues and condition numbers
        eig  = single(zeros(matricesPerFile, matrixSize));
        cond = single(zeros(matricesPerFile, matrixSize));
        
        parfor i=1:matricesPerFile

            % Generator a random matrix
            A = generator();
            
            % Compute eigenvalues and condition numbers w.r.t. eigenvalues
            [~, D, s] = condeig(A);
            
            eig(i,:)  = single(diag(D));
            cond(i,:) = single(s);
            
        end
        
        % ---------------
        % Save the data
        generatorStr = func2str(generator);
        filename = [dataDir, filenamePrefix, '_', num2str(k), '.mat'];
        save(filename, 'eig', 'cond', 'matrixSize', 'generatorStr');
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
    % writeReadMe                                                         %
    %                                                                     %
    % Write a readme file in the dataDir with information about when      %
    % the data was created, what was used to create the data, etc.        %
    % ------------------------------------------------------------------- %
    function writeReadMe()
        file = fopen([dataDir, 'README.txt'], 'w');
        fprintf(file, ['Last Update: ',datestr(now,'mmmm dd/yyyy HH:MM:SS AM'),'\n\n']);
        fprintf(file, ['Each file contains data about the eigenvalues and condition numbers w.r.t. eigenvalues for ', num2str(matricesPerFile), ' random ', num2str(matrixSize), 'x', num2str(matrixSize)]);
        fprintf(file, [' matrices where the random matrices are sampled from the ', func2str(generator), ' function handle\n']);
        fprintf(file, [filenamePrefix, '_*[0-9].mat ........ All eigenvalues and corresponding condition numbers (each file contains ', num2str(matricesPerFile*matrixSize), ' lines)\n']);
        fclose(file);
    end
    
    
    % ------------------------------------------------------------------- %
    % getStartFileIndex                                                   %
    %                                                                     %
    % Determine the index of the last file created.                       %
    %                                                                     %
    % OUTPUT                                                              %
    %   Integer, one more than the last file in the data directory with   %
    %   the filenamePrefix prefix.                                        %
    %                                                                     %
    % TO DO                                                               %
    %   Clean this up...                                                  %
    % ------------------------------------------------------------------- %
    function imax = getStartFileIndex()
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
    %       startFileIndex                                                %
    %       numFiles                                                      %
    %       matricesPerFile                                               %
    %       filenamePrefix                                                %
    %                                                                     %
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
        
        optionNames = struct('startFileIndex', 'startFileIndex', ...
                             'numFiles', 'numFiles', ...
                             'matricesPerFile', 'matricesPerFile', ...
                             'filenamePrefix', 'filenamePrefix');
        
        if ~all(ismember(fnames, fieldnames(optionNames)))
            error('processData:InvalidOption',  ...
                  'Invalid option provided');
        end
        
        % Default values
        startFileIndex  = -1;
        numFiles        = 1;
        matricesPerFile = floor(1e6/matrixSize);
        filenamePrefix  = 'BHIME';
        
        % numFiles
        if isfield(options, optionNames.numFiles)
            numFiles = options.numFiles;
            
            % Check that numFiles is a positive integer
            if ~((numFiles>0)&(mod(numFiles,1)==0))
                error('numFiles option must be a positive integer');
            end
        end
        
        % matricesPerFile
        if isfield(options, optionNames.matricesPerFile)
            matricesPerFile = options.matricesPerFile;
            
            % Check that matricesPerFile is a positive integer
            if ~((matricesPerFile>0)&(mod(matricesPerFile,1)==0))
                error('matricesPerFile option must be a positive integer');
            end
        end
        
        % filenamePrefix
        if isfield(options, optionNames.filenamePrefix)
            filenamePrefix = options.filenamePrefix;
        end
        
        % startFileIndex
        if isfield(options, optionNames.startFileIndex)
            startFileIndex = options.startFileIndex;
            
            % Check that startFileIndex is a positive integer
            if ~((startFileIndex>0)&(mod(startFileIndex,1)==0))
                error('startFileIndex option must be a positive integer');
            end
            
        else
            startFileIndex = getStartFileIndex();
        end
        
        opts = struct('startFileIndex', startFileIndex, ...
                      'numFiles', numFiles, ...
                      'matricesPerFile', matricesPerFile, ...
                      'filenamePrefix', filenamePrefix);
    end
    
end