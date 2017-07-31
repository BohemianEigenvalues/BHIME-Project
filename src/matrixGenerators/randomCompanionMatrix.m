% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Nov. 7/2016                                                 %
%                                                                         %
% This function will generate a random matrix companion matrix for a      %
% monic polynomial of degree = n+1 with coefficients sampled uniformly    %
% from the population vector.                                             %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector to sample coefficients of characteristic        %
%                  polynomial from                                        %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn companion matrix                                               %
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
function A = randomCompanionMatrix(population, n)
    
    p = [1, randsample(population, n, true)];
    A = compan(p);
    
end