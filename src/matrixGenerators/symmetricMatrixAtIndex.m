%-------------------------------------------------------------------------
% AUTHOR .... Steven E. Thornton
% EMAIL ..... sthornt7@uwo.ca
% UPDATED ... Mar. 5/2016
%
% This function generates the symmetricsquare matrix with entries from a 
% population vector such that:
%   - The population vector is of length b
%   - The max number of unique matrices is N = b^(n(n+1)/2)
%   - Let R = i mod N in base b
%   - Reverse R and pad on the left with zeros to give a number with 
%	  n(n+1)/2 digits
%   - This number then gives the indices of the population vector for each 
%     element of the output matrix
%
% INPUT
%	i ............ Index to select matrix at
%	n ............ Size of matrix
%	population ... Vector of values to be used in matrix
%
% OUTPUT
%   A symmetric square matrix where the entries are sampled from the 
%	population and determined uniquely for each value of i mod N.
%-------------------------------------------------------------------------
function A = matrixAtIndex(i, n, population)

	% Number of values in population vector
	popsize = length(population);

	% Number of different matrices belonging to this class
	% i.e. number of different nxn matrices with entries from population
	% vector
	classSize = popsize^(n*(n+1)/2);

	% Make input i value between 0 and classSize-1
	imod = mod(i, classSize);

	% Convert to base popsize to use for indexing
	idx = dec2base(imod, popsize)-'0';

	% Pad the idx vector with zeros
	idxsize = length(idx);
	idx = [zeros(1, (n*(n+1)/2) - idxsize), idx];
	idx = fliplr(idx);

	% Get population values at index and make an nxn matrix
	l = population(idx+1);
	
	A = triu(ones(n)); 
	A(A==1) = l;
	A = A + tril(A.',-1);

end