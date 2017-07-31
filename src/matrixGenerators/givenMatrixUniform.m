% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Apr. 5/2016                                                 %
%                                                                         %
% This function will take a sqaure matrix and replace specified entries   %
% by random continuous uniform random values.                             %
%                                                                         %
% INPUT                                                                   %
%   B ......... A nxn matrix                                              %
%   entries ... Matrix with 2 columns and at most n^2 rows                %
%   a ......... Left endpoint for continuous uniform sample               %
%   b ......... Right endpoint for continuous uniform sample              %
%                                                                         %
% OUTPUT                                                                  %
%   The input matrix B where the entries of B (specified by the entries   %
%   matrix) are sampled from a continuous uniform distribuion on [a,b]    %
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
function A = givenMatrixUniform(B, entries, a, b)
    
    n = size(B, 1);
    
    % Convert entries to column indexing
    idx = (entries(:,1) - 1) + n*(entries(:,2)-1) + 1;
    
    values = a + (b-a)*rand(size(entries, 1), 1);
    
    A = B;
    A(idx) = values;
    
end