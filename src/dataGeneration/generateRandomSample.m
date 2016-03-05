% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Mar. 5/2016                                                 %
%                                                                         %
% This function will generate .mat files containing the eigenvalues and   %
% their respective condition numbers for a sample of matrices. The        %
% matrices are sampled from a user provided function.                     %
%   - A readme file will be automatically generated when this function    %
%     is called                                                           %
%   - A new directory 'data/' in the input workingDir will be created     %
%   - Files will be named all_i.m where i starts at 1                     %
%                                                                         %
% INPUT                                                                   %
%   generator ......... A function handle that takes no input and returns %
%                       a random square matrix, must always return a      %
%                       matrix of the same size                           %
%   workingDir ........ The directory where files should be written       %
%                                                                         %
% OPTIONS                                                                 %
%   Options input should be a struct                                      %
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
    end
    
    % Check that generator is a function handle
    if ~isa(generator, 'function_handle')
        error('generator must be a function handle');
    end
    
    % Make data directory if it doesn't exist
    if workingDir(end) == '/'
        dataDir = [workingDir, 'Data/'];
    else
        dataDir = [workingDir, '/Data/'];
    end
    mkdir_if_not_exist(dataDir);
    
    % Convert the function handle to a string for readme
    generatorStr = func2str(generator);
    
    % Determine matrix size
    matrixSize = max(size(generator()));
    
    % Option processing
    startFileIndex = getStartFileIndex(dataDir);
    numFiles = 1;
    matricesPerFile = floor(1e6/matrixSize);
    
    if nargin == 3
        % startFileIndex
        if isfield(options,'startFileIndex')
            startFileIndex = options.startFileIndex;
        end
        
        % numFiles
        if isfield(options,'numFiles')
            numFiles = options.numFiles;
        end
        
        % matricesPerFile
        if isfield(options,'matricesPerFile')
            matricesPerFile = options.matricesPerFile;
        end
    end
    
    % Write a readme file
    writeReadMe(matrixSize, matricesPerFile, dataDir, generatorStr);
    
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
        filename = [dataDir, 'all_', num2str(k), '.mat'];
        save(filename, 'eig', 'cond', 'matrixSize', 'generatorStr');
        % ---------------
        
        % Timing
        tElapsed = toc(t2);  % timer 2
        fprintf('Loop time is %.3f seconds\n', tElapsed);
        
    end
    
    % Timing
    averageTime = toc/numFiles;
    fprintf('\nAverage loop time is %.3f seconds\n', averageTime);
    
    % Function to write a readme file
    function writeReadMe(matrixSize, matricesPerFile, dataDir, generatorStr)
        file = fopen([dataDir, 'README.txt'], 'w');
        fprintf(file, ['Last Update: ',datestr(now,'mmmm dd/yyyy HH:MM:SS AM'),'\n\n']);
        fprintf(file, ['Each file contains data about the eigenvalues and condition numbers w.r.t. eigenvalues for ', num2str(matricesPerFile), ' random ', num2str(matrixSize), 'x', num2str(matrixSize)]);
        fprintf(file, [' matrices where the random matrices are sampled from the ', generatorStr, ' function handle\n']);
        fprintf(file, ['all_*[0-9].mat ........ All eigenvalues and corresponding condition numbers (each file contains ', num2str(matricesPerFile*matrixSize), ' lines)\n']);
        fclose(file);
    end
    
    % Determine the startfile index
    % This needs to be cleaned up...
    function imax = getStartFileIndex(dataDir)
        folderInfo = dir(dataDir);
        
        imax = 0;
        
        for k=1:length(folderInfo)
            filename = folderInfo(k).name;
            
            [pathstr,name,ext] = fileparts(filename);
            
            if ext == '.mat'
                if numel(name) > 4
                    if name(1:4) == 'all_'
                        i = str2num(name(5));
                        if i > imax
                            imax = i;
                        end
                    end
                end
            end
        end
        
        % Increment by 1
        imax = imax + 1;
        
    end
    
end