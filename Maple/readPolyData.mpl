# ----------------------------------------------------------------------- #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ...  Oct. 3/2016                                                #
#                                                                         #
# Read in polynomials from a text file where each line of the file is of  #
# the form:                                                               #
#   count cn, ..., c1, c0                                                 #
# where count is the number of times the polynomial in the current line   #
# appears. c0, c1, ..., cn are the coefficients of the polynomial where   #
# the polynomial takes the form:                                          #
#   c0 + c1x + c2x^2 + ... + cn x^n           (if isMonic = false)        #
#   c0 + c1x + c2x^2 + ... + cn x^(n-1) + x^n (if isMonic = true)         #
#                                                                         #
# INPUT                                                                   #
#   filename .... Name (and location) of file.                            #
#   v ........... Variable for the output polynomials.                    #
#   isMonic ..... Optional, false by default.                             #
#                 When true, the file is assumed to only                  #
#                 contain the first n-1 coefficients as the coefficient   #
#                 of the highest degree term is 1. When false, the file   #
#                 is assumed to contain all coefficients.                 #
#   keepCount ... Optional, false by default.                             #
#                 See output.                                             #
#   delimiter ... Optional, "," by default.                               #
#                 Character that separates the entries in the file.       #
#                                                                         #
# OUTPUT                                                                  #
#   If keepCount = true:                                                  #
#       Two arrays, the first contains the polynomials and the second     #
#       contains their counts.                                            #
#   If keepCount = false:                                                 #
#       An array of polynomials                                           #
#                                                                         #
# ASSUMPTIONS                                                             #
#   - The values in the file are all assumed to be integers.              #
# ----------------------------------------------------------------------- #
readPolyData := proc(filename::string, 
                     v::symbol, 
                    {isMonic::truefalse := false, 
                     keepCount::truefalse := false, 
                     delimiter2::string := ","}, $)
    
    local data::Matrix, 
      numPolys::posint, 
             n::posint, 
   outputPolys::Array, 
         count::Array, 
             i::posint, 
             j::posint, 
             p::polynom;
    
    # Read the file into a matrix
    data := ImportMatrix(filename, 'delimiter' = delimiter2,
                                   'datatype' = integer);
    
    # Number of polynomials
    numPolys := LinearAlgebra:-RowDimension(data);
    
    # Degree of the polynomials
    if isMonic then
        n := LinearAlgebra:-ColumnDimension(data) - 1;
    else
        n := LinearAlgebra:-ColumnDimension(data) - 2;
    end if;
    
    # Array for polynomial output
    outputPolys := Array(1..numPolys, 'datatype'=polynom(integer, v));
    
    # Array for count output
    if keepCount then
        count := Array(1..numPolys, 'datatype'=nonnegint);
    end if;
    
    for i to numPolys do
        
        if isMonic then
            p := v^n;
            for j from n-1 by -1 to 0 do 
                p := p + v^j*data[i, n-j+1];
            end do;
        else
            p := 0;
            for j from n by -1 to 0 do
                p := p + v^j*data[i, n-j+2]
            end do;
        end if;
        
        outputPolys[i] := p;
        
        if keepCount then
            count[i] := data[i,1];
        end if:
        
    end do;
    
    if keepCount then
        return outputPolys, count;
    else
        return outputPolys;
    end if;
    
end proc: