% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2016)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Mar. 5/2016                                                 %
%                                                                         %
% This function will generate a random symmetric matrix with entries      %
% sampled from a list of a given output size.                             %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to be used in matrix                  %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn symmetric matrix with entries randomly sampled from the        %
%   population vector.                                                    %
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
function A = randomSymmetricMatrix(population, n)

    l = randsample(population, n*(n+1)/2, true);
    
    A = triu(ones(n)); 
    A(A==1) = l;
    A = A + tril(A.',-1);

end