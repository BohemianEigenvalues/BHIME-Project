# ----------------------------------------------------------------------- #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 26/2016                                                #                                                                         
#                                                                         #
# Get the ith matrix from a class of square matrices where the entries    #
# are from a given list                                                   #
#                                                                         #
# INPUT                                                                   #
#   i ... Index for matrix to return                                      #
#   n ... Size of matrix                                                  #
#   S ... List to select entries from                                     #
#                                                                         #
# OUTPUT                                                                  #
#   An nxn matrix with entries in S.                                      #
# ----------------------------------------------------------------------- #
getMatrixAtIndex := proc(i, n::posint, S::list, $)::Matrix;

    local matList, c;

    # Build the matrix
    matList := NULL:

    c := i:
    to n^2 do
        matList := matList, S[irem(c, nops(S))+1];
        c := iquo(c, nops(S));
    end do;

    return Matrix(n, n, [matList]);

end proc: