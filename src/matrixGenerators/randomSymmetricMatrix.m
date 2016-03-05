%-------------------------------------------------------------------------
% AUTHOR .... Steven E. Thornton
% EMAIL ..... sthornt7@uwo.ca
% UPDATED ... Mar. 5/2016
%
% This function will generate a random symmetric matrix with entries
% sampled from a list of a given output size.
%
% INPUT
%	population ... Vector of values to be used in matrix
%	n ............ Size of matrix
%
% OUTPUT
%   An nxn symmetric matrix with entries randomly sampled from the 
%   population vector
%-------------------------------------------------------------------------
function A = randomSymmetricMatrix(population, n)

    l = randsample(population, n*(n+1)/2, true);
    
    A = triu(ones(n)); 
    A(A==1) = l;
    A = A + tril(A.',-1);

end