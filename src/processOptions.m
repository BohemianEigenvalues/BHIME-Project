% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Apr. 5/2016                                                 %
%                                                                         %
% Process the options input options struct. If an option is not in the    %
% options struct the default value is used.                               %
%                                                                         %
% INPUT                                                                   %
%   options ... (struct) contains keys corresponding to the options       %
%                                                                         %
% OUTPUT                                                                  %
%   A struct opts with keys                                               %
%       startFileIndex .... posint                                        %
%       numDataFiles ...... posint                                        %
%       numProcessFiles ... posint                                        %
%       matricesPerFile ... posint                                        %
%       filenamePrefix .... str, valid filename                           %
%       height ............ posint                                        %
%       margin ............ struct(left, right, bottom, top), double      %
%       outputFileType .... {'txt', 'mat'}                                %
%       symmetry .......... bool                                          %
%       map ............... Function handle, one input value.             %
%                           double -> double                              %
%       backgroundColor ... [double, double, double] in [0, 1]            %
%                                                                         %
% TO DO                                                                   %
%   - Add type checking for options                                       %
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
function opts = processOptions(options)
    
    % Check the number of input arguments
    if nargin ~= 1
        error('processOptions:WrongNumberofArgs', ...
              'processOptions expects one argument.');
    end
    
    % Check that options is a struct
    if ~isstruct(options)
        error('processOptions:InvalidOptionsStruct', ...
              'options argument must be a structured array');
    end
    
    optNames = struct( 'startFileIndex', 'startFileIndex', ...
                         'numDataFiles', 'numDataFiles', ...
                      'numProcessFiles', 'numProcessFiles', ...
                      'matricesPerFile', 'matricesPerFile', ...
                       'filenamePrefix', 'filenamePrefix', ...
                               'height', 'height', ...
                               'margin', 'margin', ...
                       'outputFileType', 'outputFileType', ...
                             'symmetry', 'symmetry', ...
                                  'map', 'map', ...
                      'backgroundColor', 'backgroundColor');
    
    fnames = fieldnames(options);
    
    if ~all(ismember(fnames, fieldnames(optNames)))
        error('processOptions:InvalidOption',  ...
              'Invalid option provided');
    end
    
    opts = struct();
    
    % Default values
    opts.startFileIndexIsSet  = false;
    opts.startFileIndex       = 1;
    opts.numDataFiles         = 1;
    opts.numProcessFilesIsSet = false;
    opts.numProcessFiles      = 1;
    opts.matricesPerFileIsSet = false;
    opts.matricesPerFile      = 1e6;
    opts.filenamePrefix       = 'BHIME';
    opts.height               = 1001;
    opts.marginIsSet          = false;
    opts.margin               = struct();
    opts.outputFileType       = 'mat';
    opts.symmetry             = false;
    opts.map                  = @(z) z;
    opts.backgroundColor      = [0, 0, 0];
    
    % startFileIndex -------------------
    if isfield(options, optNames.startFileIndex)
        
        opts.startFileIndex = options.startFileIndex;
        opts.startFileIndexIsSet = true;
        
        % Check that startFileIndex is a positive integer
        if ~((startFileIndex>0)&(mod(startFileIndex,1)==0))
            error('startFileIndex option must be a positive integer');
        end
    end
    
    % numDataFiles ---------------------
    if isfield(options, optNames.numDataFiles)
        
        opts.numDataFiles = options.numDataFiles;
        
        % Check that numFiles is a positive integer
        if ~((opts.numDataFiles>0)&(mod(opts.numDataFiles,1)==0))
            error('numDataFiles option must be a positive integer');
        end
        
    end
    
    % numProcessFiles ------------------
    if isfield(options, optNames.numProcessFiles)
        
        opts.numProcessFiles = options.numProcessFiles;
        opts.numProcessFilesIsSet = true;
        
        % Check that numFiles is a positive integer
        if ~((opts.numProcessFilesIsSet>0)&(mod(opts.numProcessFilesIsSet,1)==0))
            error('numProcessFilesIsSet option must be a positive integer');
        end
        
    end
    
    % matricesPerFile ------------------
    if isfield(options, optNames.matricesPerFile)
        
        opts.matricesPerFile = options.matricesPerFile;
        opts.matricesPerFileIsSet = true;
        
        % Check that matricesPerFile is a positive integer
        if ~((opts.matricesPerFile>0)&(mod(opts.matricesPerFile,1)==0))
            error('matricesPerFile option must be a positive integer');
        end
        
    end
    
    % filenamePrefix -------------------
    if isfield(options, optNames.filenamePrefix)
        opts.filenamePrefix = options.filenamePrefix;
        
        % TO DO: Add type checking (string, valid filename prefix)
        
    end
    
    % height ---------------------------
    if isfield(options, optNames.height)
        opts.height = options.height;
        
        % TO DO: Add type checking (positive integer)
        
    end
    
    % margin ---------------------------
    if isfield(options, optNames.margin)
        
        opts.margin = options.margin;
        opts.marginIsSet = true;
        
        % TO DO: Add type checking (struct: left, right, bottom, top)
    end
    
    % outputFileType -------------------
    if isfield(options, optNames.outputFileType)
        
        opts.outputFileType = options.outputFileType;
        
        % TO DO: Add type checking (txt or mat)
        
    end
    
    % symmetry -------------------------
    if isfield(options, optNames.symmetry)
        
        opts.symmetry = options.symmetry;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % map ------------------------------
    if isfield(options, optNames.map)
        
        opts.map = options.map;
        
        % TO DO: Add type checking (function handle with 1 input)
        
    end
    
    % backgroundColor ------------------
    if isfield(options, optNames.backgroundColor)
        
        opts.backgroundColor = options.backgroundColor;
        
        % TO DO: Add type checking (vector of 3 values in [0,1])
        
    end
    
end