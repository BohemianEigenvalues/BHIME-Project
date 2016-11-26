# ----------------------------------------------------------------------- #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 3/2016                                                 #
#                                                                         #
# Count the number of derogatory matrices (charPoly <> minPoly) in a      #
# class of matrices.                                                      #
#                                                                         #
# INPUT                                                                   #
#   n ... Matrix size                                                     #
#   S ... List of set of entries                                          #
#                                                                         #
# OUTPUT                                                                  #
#   If output = "table":                                                  #
#       A table is returned where the keys are the characteristic         #
#       polynomials and the values are their densities.                   #
#   If output = "set":                                                    #
#       A set of the characteristic polynomials is returned.              #
# ----------------------------------------------------------------------- #
countDerogarory := proc(n::posint, S::{list, set}, $)::nonnegint;
    
    local classSize, count, i, A, charPoly, minPoly;
    
    # Number of matrices in the given class
    classSize := nops(S)^(n^2);
    
    count := 0;
    
    for i to classSize do
    
        A := getMatrixAtIndex(i, n, S):
        charPoly := expand(LinearAlgebra:-CharacteristicPolynomial(A, 'x'));
        minPoly := expand(LinearAlgebra:-MinimalPolynomial(A, 'x'));
        
        if abs(normal(charPoly/minPoly)) <> 1 then
            count := count + 1;
        end if;
        
    end do;
    
    return count;
    
end proc:
