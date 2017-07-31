% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Jul. 31/2017                                                %
%                                                                         %
% This function will generate a random non-normal upper-hessenberg matrix %
% with a Toeplitz structure and a zero main diagonal where the entries    %
% are sampled from a given list.                                          %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to sample entries from                %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn non-normal upper-hessenberg Toeplitz matrix with zero main     %
%   diagonal with entries sampled from the population vector.             %
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
function A = randomNonNormalUpperHessenbergToeplitzZeroDiagMatrix(population, n)
    while true
        % Get a Toeplitz matrix
        A = randomToeplitzMatrix(population, n);
        
        % Make it upper-hessenberg
        A = triu(A,-1);
        
        % Set the main diagonal to zero
        A(eye(size(A))~=0)=0;
        
        % Check if it is non-normal
        ANorm = A*(A') - (A')*A;
        if ~all(ANorm(:) == 0)
            return
        end
    end
end