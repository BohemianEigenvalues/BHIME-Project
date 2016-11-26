# ----------------------------------------------------------------------- #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 3/2016                                                 #
#                                                                         #
# Compute the square-free decomposition of a polynomial and return the    #
# product of the square-free terms.                                       #
#                                                                         #
# INPUT                                                                   #
#   p ... Polynomial with integer coefficients                            #
#                                                                         #
# OUTPUT                                                                  #
#   A polynomial with no multiple roots.                                  #
# ----------------------------------------------------------------------- #
makePolySquareFree := proc(p::polynom, $)::polynom;
    
    return mul(map(proc (a) options operator, arrow; a[1] end proc, sqrfree(p)[2]));
    
end proc: