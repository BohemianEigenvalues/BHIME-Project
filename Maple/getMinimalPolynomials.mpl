# ----------------------------------------------------------------------- #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 3/2016                                                 #
#                                                                         #
# This function will compute the set of all minimal polynomials for       #
# matrices where the entries come from a given set.                       #
#                                                                         #
# INPUT                                                                   #
#   n ... Matrix size                                                     #
#   S ... List of set of entries                                          #
#                                                                         #
# OPTIONS                                                                 #
#   printEvery ... Default = 10,000                                       #
#                  How often to print a message                           #
#   output ....... Default = "table"                                      #
#                  Either "table" or "set"                                #
#                                                                         #
# OUTPUT                                                                  #
#   If output = "table":                                                  #
#       A table is returned where the keys are the minimal polynomials    #
#       and the values are their densities.                               #
#   If output = "set":                                                    #
#       A set of the minimal polynomials is returned.                     #
# ----------------------------------------------------------------------- #
getMinimalPolynomials := proc(n::posint, 
                              S::{list, set}, 
                             {printEvery::posint := 10000,
                              output::string := "table"}, 
                              $)::{table, set};

    local classSize, T, i, A, p;
    
    # Number of matrices in the given class
    classSize := nops(S)^(n^2);

    if output = "table" then

        T := table();

        for i to classSize do

            if modp(i, printEvery) = 0 then
                printf("%d of %d\n", i, classSize);
            end if;

            A := getMatrixAtIndex(i, n, S):
            p := expand(LinearAlgebra:-MinimalPolynomial(A, 'x'));

            if p in {op(ListTools:-Flatten([indices(T)]))} then
                T[p] := T[p] + 1;
            else
                T[p] := 1;
            end if;

        end do:

    elif output = "set" then

        T := {};

        for i to classSize do

            if modp(i, printEvery) = 0 then
                printf("%d of %d\n", i, classSize);
            end if;

            A := getMatrixAtIndex(i, n, S):
            p := expand(LinearAlgebra:-MinimalPolynomial(A, 'x'));

            T := `union`(T, {p});

        end do:

    else
        error("Incorrect output option type. Expected Table or Set.");
    end if;

    return T;

end proc:

