% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Jul. 31/2017                                                %
%                                                                         %
% This function will generate a random Hamiltonian matrix where the       %
% entries are sampled from a given list.                                  %
%                                                                         %
% INPUT                                                                   %
%   population ... Vector of values to sample entries from                %
%   n ............ Size of matrix                                         %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn Hamiltonian matrix where the entries are randomly sampled from %
%   the population vector. i.e. A matrix of the form:                     %
%           [A B                                                          %
%            C -A'];                                                      %
%   where B and C are symmetric matrices. A, B, and C are all of square   %
%   matrices of the same size.                                            %
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
function H = randomHamiltonianMatrix(population, n)
    
    A = randsample(population, (n/2)^2, true);
    A = reshape(A, [n/2, n/2]);
    
    l = randsample(population, (n/2)*((n/2) + 1)/2, true);
    B = triu(ones(n/2)); 
    B(B==1) = l;
    B = B + tril(B.',-1);
    
    l = randsample(population, (n/2)*((n/2) + 1)/2, true);
    C = triu(ones(n/2)); 
    C(C==1) = l;
    C = C + tril(C.',-1);
   
    H = [A, B;
         C, -A'];
    
end