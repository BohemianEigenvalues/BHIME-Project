% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Apr. 22/2016                                                %
%                                                                         %
% This function will generate a random Toeplitz matrix of a given size    %
% where the entries are sampled uniformly from a given list.              %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to be used in matrix                  %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn Toeplitz matrix with entries randomly sampled from the         %
%   population vector                                                     %
%                                                                         %
% REFERENCES                                                              %
%   wikipedia.org/wiki/Toeplitz_matrix                                    %
%   mathworld.wolfram.com/ToeplitzMatrix.html                             %
%   mathworks.com/help/matlab/ref/toeplitz.html                           %
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
function A = randomToeplitzMatrix(population, n)

    vec = randsample(population, 2*n-1, true);
    
    A = toeplitz(vec(n:-1:1), vec(n:end));

end