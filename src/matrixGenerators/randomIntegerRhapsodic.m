% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Jul. 31/2017                                                %
%                                                                         %
% This function will generate a random rhapsodic matrix with entries      %
% sampled from a given list of integers. A rhapsodic matrix is a matrix   %
% whose inverse belongs to the same family of Bohemian matrices that the  %
% original matrix was sampled from.                                       %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to sample entries from. Must only     %
%                  contain integers.                                      %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn rhapsodic matrix.                                              %
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
function A = randomIntegersRhapsodic(population, n)
    while true
        A = randsample(population, n^2, true);
        A = reshape(A, [n, n]);
        if abs(det(A)) == 1
            Ainv = inv(A);
            if isempty(setdiff(Ainv, population))
                return
            end
        end
    end
end