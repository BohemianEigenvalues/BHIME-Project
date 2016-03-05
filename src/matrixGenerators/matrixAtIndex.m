% ------------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                       %
% EMAIL ..... sthornt7@uwo.ca                                               %
% UPDATED ... Mar. 5/2016                                                   %
%                                                                           %
% This function generates the square matrix with entries from a population  %
% vector such that:                                                         %
%   - The population vector is of length b                                  %
%   - The max number of unique matrices is N = b^(n^2)                      %
%   - Let R = i mod N in base b                                             %
%   - Reverse R and pad on the left with zeros to give a number with n^2    %
%     digits                                                                %
%   - This number then gives the indices of the population vector for each  %
%     element of the output matrix                                          %
%                                                                           %
% INPUT                                                                     %
%	i ............ Index to select matrix at                                %
%	n ............ Size of matrix                                           %
%	population ... Vector of values to be used in matrix                    %
%                                                                           %
% OUTPUT                                                                    %
%   A square matrix where the entries are sampled from the population and   %
%   determined uniquely for each value of i mod N.                          %
%                                                                           %
% LICENSE                                                                   %
%	This program is free software: you can redistribute it and/or modify    %
%	it under the terms of the GNU General Public License as published by    %
% 	the Free Software Foundation, either version 3 of the License, or       %
% 	any later version.                                                      %
%                                                                           %
% 	This program is distributed in the hope that it will be useful,         %
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           %
% 	GNU General Public License for more details.                            %
%                                                                           %
% 	You should have received a copy of the GNU General Public License       %
% 	along with this program.  If not, see http://www.gnu.org/licenses/.     %
% ------------------------------------------------------------------------- %
function A = matrixAtIndex(i, n, population)
    
    % Number of values in population vector
    popsize = length(population);

    % Number of different matrices belonging to this class
    % i.e. number of different nxn matrices with entries from population
    % vector
    classSize = popsize^(n^2);
    
    % Make input i value between 0 and classSize-1
    imod = mod(i, classSize);
    
    % Convert to base popsize to use for indexing
    idx = dec2base(imod, popsize)-'0';
    
    % Pad the idx vector with zeros
    idxsize = length(idx);
    idx = [zeros(1, n^2 - idxsize), idx];
    idx = fliplr(idx);
    
    % Get population values at index and make an nxn matrix
    A = reshape(population(idx+1), [n,n])';
    
end