% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Mar. 11/2016                                                %
%                                                                         %
% This function will generate .mat files containing the eigenvalues and   %
% their respective condition numbers for a all matrices matrices from the %
% input generating function.                                              %
%   - A readme file will be automatically generated when this function    %
%     is called                                                           %
%   - A new directory 'Data' in the input workingDir will be created      %
%   - Files will be named BHIME_i.mat where i starts at 1 (you can        %
%     customize the prefix by using the filenamePrefix option)            %
%                                                                         %
% INPUT                                                                   %
%   workingDir ........ The directory where files should be written       %
%   generator ......... A function handle that takes a positive integer   %
%                       as input and returns a square matrix that is      %
%                       unique for that input integer. Must always return %
%                       a matrix of the same size.                        %
%   numMatrices ....... The total number of unique matrices the generator %
%                       can generate.                                     %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
%   filenamePrefix .... Default = BHIME                                   %
%                       The prefix for the .mat files that contain the    %
%                       eigenvalues and their condition numbers           %
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
function generateAllMatrices(workingDir, generator, numMatrices, options)
    
    % Check number of input arguments
    if nargin > 4
        error('generateRandomSample:TooManyInputs', ...
              'requires at most 4 input arguments');
    elseif nargin < 3
        error('generateRandomSample:NotEnoughInputs', ...
              'requires at least 3 input arguments');
    elseif nargin == 3
        options = struct()
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
    matrixSize = max(size(generator(1)));
    
    % Process the options
    opts = processOptions();
    matricesPerFile = opts.matricesPerFile;
    filenamePrefix = opts.filenamePrefix;
    
    % Write a readme file
    writeReadMe();
    
    % Set random seed based on the current time 
    s = RandStream('mt19937ar', 'Seed', mod(int64(now*1e6), 2^32));
    RandStream.setGlobalStream(s);
    
    % Master timer
    tic;
    
    % Split into groups of size matricesPerFile
    numOuterLoops = ceil(numMatrices/matricesPerFile);
    
    for i=1:numOuterLoops
        
        % Timer 2
        t2 = tic;
        
        fprintf('\n%d of %d\n',i,numOuterLoops);
        
        % Determine starting and ending 
        jStart = (i-1)*matricesPerFile + 1;
        jEnd = i*matricesPerFile;
        
        if jEnd > numMatrices
            % Check if we're on the last file
            jEnd = numMatrices;
            matricesPerFile = jEnd - jStart + 1;
        end
        
        % Vectors to store eigenvalues and condition numbers
        eig  = single(zeros(matricesPerFile, matrixSize));
        cond = single(zeros(matricesPerFile, matrixSize));
        
        parfor j=1:matricesPerFile
        
            % Offset counter
            jj = jStart + j - 1;
            
            % Get the matrix
            A = generator(jj);
            
            % Compute eigenvalues and condition numbers w.r.t. eigenvalues
            [~, D, s] = condeig(A);
            
           eig(j,:)  = single(diag(D));
           cond(j,:) = single(s);
        
        end
        
        % ---------------
        % Save the data
        generatorStr = func2str(generator);
        filename = [dataDir, filenamePrefix, '_', num2str(i), '.mat'];
        save(filename, 'eig', 'cond', 'matrixSize', 'generatorStr');
        % ---------------
        
        % Timing
        tElapsed = toc(t2);  % timer 2
        fprintf('Loop time is %.3f seconds\n', tElapsed);
        
    end
    
    % Timing
    averageTime = toc/numOuterLoops;
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
        t = datestr(now, 'mmmm dd/yyyy HH:MM:SS AM');
        fprintf(file, ['Last Update: ', t, '\n\n']);
        fprintf(file, 'Each file contains data about the (right) ');
        fprintf(file, 'eignenvalues and condition numbers w.r.t. the ');
        fprintf(file, 'eigenvalues for ');
        fprintf(file, [num2str(matricesPerFile), ' ']);
        fprintf(file, [num2str(matrixSize), 'x', num2str(matrixSize)]);
        fprintf(file, ' matrices where the matrices are generated from ');
        fprintf(file, ['the ' func2str(generator), ' function handle\n']);
        fprintf(file, [filenamePrefix, '_*[0-9].mat ... All eigenvalues']);
        fprintf(file, ' and corresponding condition numbers (each file ');
        fprintf(file, ['contains ' num2str(matricesPerFile*matrixSize)]);
        fprintf(file, [' lines\n']);
        fclose(file);
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
        
        optionNames = struct('matricesPerFile', 'matricesPerFile', ...
                             'filenamePrefix', 'filenamePrefix');
        
        if ~all(ismember(fnames, fieldnames(optionNames)))
            error('processData:InvalidOption',  ...
                  'Invalid option provided');
        end
        
        % Default values
        matricesPerFile = floor(1e6/matrixSize);
        filenamePrefix  = 'BHIME';
        
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
        
        opts = struct('matricesPerFile', matricesPerFile, ...
                      'filenamePrefix', filenamePrefix);
    end
    
    
end