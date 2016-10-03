% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... May 22/2016                                                 %
%                                                                         %
% Process the options input options struct. If an option is not in the    %
% options struct the default value is used.                               %
%                                                                         %
% INPUT                                                                   %
%   options ... (struct) contains keys corresponding to the options       %
%                                                                         %
% OUTPUT                                                                  %
%   A struct opts with keys                                               %
%       backgroundColor ............ [double, double, double] all in [0,1]%
%       colorByCond ................ bool                                 %
%       computeCond ................ bool                                 %
%       dataPrecision .............. str in {'double', 'single'}          %
%       filenamePrefix ............. str, valid filename                  %
%       height ..................... posint                               %
%       ignoreReal ................. bool                                 %
%       ignoreRealData ............. bool                                 %
%       ignoreRealTol .............. bool                                 %
%       map ........................ Function handle:                     %
%                                       (double) -> double                %
%       margin ..................... struct of doubles of the form:       %
%                                       'bottom': double                  %
%                                         'top' : double                  %
%                                        'left' : double                  %
%                                       'right' : double                  %
%       matricesPerFile ............ posint                               %
%       maxDensity ................. posint                               %
%       minDensity ................. nonnegint                            %
%       numCharPolyFiles ........... posint                               %
%       numDataFiles ............... posint                               %
%       numProcessFiles ............ posint                               %
%       outputFileType ............. str in {'txt', 'mat'}                %
%       overrideDataDir ............ str, valid directory                 %
%       overrideImagesDir .......... str, valid directory                 %
%       overrideProcessedDataDir ... str, valid directory                 %
%       startFileIndex ............. posint                               %
%       storeDataWithSymmetry ...... bool                                 %
%       symmetryIm ................. bool                                 %
%       symmetryRe ................. bool                                 %
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
    
    optNames = struct('backgroundColor', 'backgroundColor', ...
                          'colorByCond', 'colorByCond', ...
                          'computeCond', 'computeCond', ...
                        'dataPrecision', 'dataPrecision', ...
                       'filenamePrefix', 'filenamePrefix', ...
                               'height', 'height', ...
                           'ignoreReal', 'ignoreReal', ...
                       'ignoreRealData', 'ignoreRealData', ...
                        'ignoreRealTol', 'ignoreRealTol', ...
                                  'map', 'map', ...
                               'margin', 'margin', ...
                      'matricesPerFile', 'matricesPerFile', ...
                           'maxDensity', 'maxDensity', ...
                           'minDensity', 'minDensity', ...
                     'numCharPolyFiles', 'numCharPolyFiles', ...
                         'numDataFiles', 'numDataFiles', ...
                      'numProcessFiles', 'numProcessFiles', ...
                       'outputFileType', 'outputFileType', ...
                      'overrideDataDir', 'overrideDataDir', ...
                    'overrideImagesDir', 'overrideImagesDir', ...
             'overrideProcessedDataDir', 'overrideProcessedDataDir', ...
                       'startFileIndex', 'startFileIndex', ...
                'storeDataWithSymmetry', 'storeDataWithSymmetry', ...
                           'symmetryIm', 'symmetryIm', ...
                           'symmetryRe', 'symmetryRe');
    
    fnames = fieldnames(options);
    
    if ~all(ismember(fnames, fieldnames(optNames)))
        error('processOptions:InvalidOption',  ...
              'Invalid option provided');
    end
    
    opts = struct();
    
    % Default values
    opts.backgroundColor               = [0, 0, 0];
    opts.colorByCond                   = false;
    opts.computeCond                   = false;
    opts.dataPrecision                 = 'single';
    opts.filenamePrefix                = 'BHIME';
    opts.height                        = 1001;
    opts.ignoreReal                    = false;
    opts.ignoreRealData                = false;
    opts.ignoreRealTol                 = 1e-10;
    opts.map                           = @(z) z;
    opts.margin                        = struct();
    opts.marginIsSet                   = false;
    opts.matricesPerFile               = 1e6;
    opts.matricesPerFileIsSet          = false;
    opts.maxDensity                    = 0;
    opts.maxDensityIsSet               = false;
    opts.minDensity                    = 0;
    opts.minDensityIsSet               = false;
    opts.numCharPolyFiles              = 1;
    opts.numCharPolyFilesIsSet         = false;
    opts.numDataFiles                  = 1;
    opts.numProcessFilesIsSet          = false;
    opts.numProcessFiles               = 1;
    opts.outputFileType                = 'mat';
    opts.overrideImagesDir             = '';
    opts.overrideImagesDirIsSet        = false;
    opts.overrideDataDir               = '';
    opts.overrideDataDirIsSet          = false;
    opts.overrideProcessedDataDir      = '';
    opts.overrideProcessedDataDirIsSet = false;
    opts.startFileIndexIsSet           = false;
    opts.startFileIndex                = 1;
    opts.storeDataWithSymmetry         = false;
    opts.symmetryIm                    = false;
    opts.symmetryRe                    = false;
    
    
    % backgroundColor ------------------
    if isfield(options, optNames.backgroundColor)
        
        opts.backgroundColor = options.backgroundColor;
        
        % TO DO: Add type checking (vector of 3 values in [0,1])
        
    end
    
    % colorByCond ----------------------
    if isfield(options, optNames.colorByCond)
        
        opts.colorByCond = options.colorByCond;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % computeCond ----------------------
    if isfield(options, optNames.computeCond)
        
        opts.computeCond = options.computeCond;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % dataPrecision --------------------
    if isfield(options, optNames.dataPrecision)
        
        opts.dataPrecision = options.dataPrecision;
        
        % TO DO: Add type checking {'single', 'double'}
        
    end
    
    % filenamePrefix -------------------
    if isfield(options, optNames.filenamePrefix)
        opts.filenamePrefix = options.filenamePrefix;
        
        % TO DO: Add type checking (string, valid filename prefix)
        if ~isValidString(opts.filenamePrefix)
            error('filenamePrefix option must be a string');
        end
    end
    
    % height ---------------------------
    if isfield(options, optNames.height)
        
        opts.height = options.height;
        
        % TO DO: Add type checking (positive integer)
        if ~isposint(opts.height)
            error('height option must be a positive integer');
        end
        
    end
    
    % ignoreReal -----------------------
    if isfield(options, optNames.ignoreReal)
        
        opts.ignoreReal = options.ignoreReal;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % ignoreRealData -------------------
    if isfield(options, optNames.ignoreRealData)
        
        opts.ignoreRealData = options.ignoreRealData;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % ignoreRealTol --------------------
    if isfield(options, optNames.ignoreRealTol)
        
        opts.ignoreRealTol = options.ignoreRealTol;
        
        % TO DO: Add type checking (double between 0 and 1)
        
    end
    
    % map ------------------------------
    if isfield(options, optNames.map)
        
        opts.map = options.map;
        
        % TO DO: Add type checking (function handle with 1 input)
        
    end
    
    % margin ---------------------------
    if isfield(options, optNames.margin)
        
        opts.margin = options.margin;
        opts.marginIsSet = true;
        
        % TO DO: Add type checking (struct: left, right, bottom, top)
    end
    
    % matricesPerFile ------------------
    if isfield(options, optNames.matricesPerFile)
        
        opts.matricesPerFile = options.matricesPerFile;
        opts.matricesPerFileIsSet = true;
        
        % Check that matricesPerFile is a positive integer
        if ~isposint(opts.matricesPerFile)
            error('matricesPerFile option must be a positive integer');
        end
    end
    
    % maxDensity -----------------------
    if isfield(options, optNames.maxDensity)
        
        opts.maxDensity = options.maxDensity;
        opts.maxDensityIsSet = true;
        
        % Check that matricesPerFile is a positive integer
        if ~isposint(opts.matricesPerFile)
            error('maxDensity option must be a positive integer');
        end
        
    end
    
    % minDensity -----------------------
    if isfield(options, optNames.minDensity)
        
        opts.minDensity = options.minDensity;
        opts.minDensityIsSet = true;
        
    end
    
    % numCharPolyFiles -----------------
    if isfield(options, optNames.numCharPolyFiles)
        
        opts.numCharPolyFiles = options.numCharPolyFiles;
        opts.numCharPolyFilesIsSet = true;
        
        % Check that numFiles is a positive integer
        if ~isposint(opts.numCharPolyFiles)
            error('numCharPolyFiles option must be a positive integer');
        end
        
    end
    
    % numDataFiles ---------------------
    if isfield(options, optNames.numDataFiles)
        
        opts.numDataFiles = options.numDataFiles;
        
        % Check that numFiles is a positive integer
        if ~isposint(opts.numDataFiles)
            error('numDataFiles option must be a positive integer');
        end
        
    end
    
    % numProcessFiles ------------------
    if isfield(options, optNames.numProcessFiles)
        
        opts.numProcessFiles = options.numProcessFiles;
        opts.numProcessFilesIsSet = true;
        
        % Check that numFiles is a positive integer
        if ~isposint(opts.numProcessFiles)
            error('numProcessFiles option must be a positive integer');
        end
    end
    
    % outputFileType -------------------
    if isfield(options, optNames.outputFileType)
        
        opts.outputFileType = options.outputFileType;
        
        % TO DO: Add type checking (txt or mat)
        if opts.outputFileType ~= 'txt' || opts.outputFileType ~= 'mat'
            error('outputFileType option must be either ''txt'' or ''mat''');
        end
        
    end
    
    % overrideDataDir ------------------
    if isfield(options, optNames.overrideDataDir)
        
        opts.overrideDataDir = options.overrideDataDir;
        opts.overrideDataDirIsSet = true;
        
        % TO DO: Add type checking, str, valid directory
        
    end
    
    % overrideImagesDir ----------------
    if isfield(options, optNames.overrideImagesDir)
        
        opts.overrideImagesDir = options.overrideImagesDir;
        opts.overrideImagesDirIsSet = true;
        
        % TO DO: Add type checking, str, valid directory
        
    end
    
    % overrideProcessedDataDir ---------
    if isfield(options, optNames.overrideProcessedDataDir)
        
        opts.overrideProcessedDataDir = options.overrideProcessedDataDir;
        opts.overrideProcessedDataDirIsSet = true;
        
        % TO DO: Add type checking, str, valid directory
        
    end
    
    % startFileIndex -------------------
    if isfield(options, optNames.startFileIndex)
        
        opts.startFileIndex = options.startFileIndex;
        opts.startFileIndexIsSet = true;
        
        % Check that startFileIndex is a positive integer
        if ~isposint(opts.startFileIndex)
            error('startFileIndex option must be a positive integer');
        end
    end
    
    % storeDataWithSymmetry ------------
    if isfield(options, optNames.storeDataWithSymmetry)
        
        opts.storeDataWithSymmetry = options.storeDataWithSymmetry;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % symmetryIm -----------------------
    if isfield(options, optNames.symmetryIm)
        
        opts.symmetryIm = options.symmetryIm;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % symmetryRe -----------------------
    if isfield(options, optNames.symmetryRe)
        
        opts.symmetryRe = options.symmetryRe;
        
        % TO DO: Add type checking (bool)
        
    end
    
    % ---------------------------------------------------------------------
    % Type checking functions
    % ---------------------------------------------------------------------
    
    function b = isposint(x)
        
        b = logical(0);
        
        try
            b = logical((x > 0) & (~mod(x, 1)));
        end
    end
    
    % Source: http://www.mathworks.com/matlabcentral/fileexchange/8847-dstv--datashape---type-verification/content/dstv/isstring.m
    function b = isValidString(S)
        % ISSTRING True for 1-D char arrays.
        %   ISSTRING(S) returns 1 if S is a string. A string is a 1-D char
        %   array and is therefore the char equivalent of a vector.
        %
        %   See also ISVECTOR, ARGCHK.
        %
        % @author   Ingo Lhken, ingo.loehken@gmx.de
        % @rev      DSTV ver 0.1, 25/10/2005
        b = logical(0);
        
        try
            b = logical(ischar(S) && ndims(S) == 2 && any(size(S) <= 1));
        end
    end
end