% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Jul. 31/2017                                                %
%                                                                         %
% This function will generate a random pentadiagonal matrix where the     %
% entries are sampled from a given list.                                  %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to sample entries from                %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn pentadiagonal matrix where the entries are randomly sampled    %
%   from the population vector.                                           %
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
function A = randomPentadiagonalMatrix(population, n)
    
    v1 = randsample(population, n, true);
    v2 = randsample(population, n-1, true);
    v3 = randsample(population, n-1, true);
    v4 = randsample(population, n-2, true);
    v5 = randsample(population, n-2, true);
    
    A = diag(v1) + diag(v2, 1) + diag(v3, -1) + diag(v4, 2) + diag(v5, -2);

end