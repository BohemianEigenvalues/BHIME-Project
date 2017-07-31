% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Jul. 31/2017                                                %
%                                                                         %
% This function will generate a random "inverted checkerboard" matrix     %
% where the entries are sampled from a given list.                        %
%                                                                         %
% An "inverted checkerboard" matrix contains zeros on the odd sub- and    %
% superdiagonals (i.e. the 1st, 3rd, 5th, etc. sub- and superdiagonals).  %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to sample entries from                %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn "inverted checkerboard" matrix where the entries are randomly  %
%   sampled from the population vector.                                   %
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
function A = randomInvertedCheckerboardMatrix(population, n)
    
    zeroband = zeros(n,1);
    for i=1:2:n
        zeroband(i) = 1;
    end
    
    Aidx = toeplitz(zeroband, zeroband);
    
    A = randsample(population, n^2, true);
    A = reshape(A, [n, n]);
    A(Aidx == 1) = 0;

end