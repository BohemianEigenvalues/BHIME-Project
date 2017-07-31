% ----------------------------------------------------------------------- %
% AUTHOR .... Steven E. Thornton (Copyright (c) 2017)                     %
% EMAIL ..... sthornt7@uwo.ca                                             %
% UPDATED ... Jul. 31/2017                                                %
%                                                                         %
% This function will generate a random matrix where the entries are       %
% sampled from a shifted and scaled beta distribution. Entries iid.       %
%                                                                         %
% INPUT                                                                   %
%   alpha ... Shape parameter for beta distribution                       %
%   beta .... Shape parameter for beta distribution                       %
%   shift ... Amount to shift values by                                   %
%   scale ... Amount to scale values by                                   %
%   n ....... Size of matrix                                              %
%                                                                         %
% OUTPUT                                                                  %
%   An nxn matrix where the entries are of the form:                      %
%       A_{i,j} = shift + scale*X_{i,j}                                   %
%   where X_{i,j} is sampled from a beta distribution with parameters     %
%   alpha and beta. Note that a beta distribution will return a values    %
%   between 0 and 1.                                                      %
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
function A = randomBeta(alpha, beta, shift, scale, n)
    A = shift + scale*random('beta', alpha, beta, [n,n]);
end