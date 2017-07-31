% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Jul. 31/2017                                                %
%                                                                         %
% This function will generate a random Cauchy matrix where the entries    %
% generated from two disctinct vectors of values sampled from a given     %
% list.                                                                   %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to sample vectors from                %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn Cauchy matrix where the two distinct vectors of length n are   %
%   sampled from the imput population vector, and a Cauchy matrix is      %
%   generated from these vectors.                                         %
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
function A = randomCauchyMatrix(population, n)
    
    x = randsample(population, n, false);
    
    newpop = setdiff(population, x);
    
    y = randsample(newpop, n, false);
    
    A = bsxfun(@(a, b) 1./(a-b), x, y');
    
end