%-------------------------------------------------------------------------
% AUTHOR .... Steven E. Thornton
% EMAIL ..... sthornt7@uwo.ca
% UPDATED ... Mar. 5/2016
%
% This function will generate a random matrix with entries
% sampled from a list of a given output size.
%
% INPUT
%	population ... Vector of values to be used in matrix
%	n ............ Size of matrix
%
% OUTPUT
%   An nxn matrix with entries randomly sampled from the population vector
%-------------------------------------------------------------------------
function A = randomMatrix(population, n)

	A = randsample(population, n^2, true);
	A = reshape(A, [n, n]);

end