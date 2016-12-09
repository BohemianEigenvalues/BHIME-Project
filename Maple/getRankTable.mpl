getRankTable := proc(n, S)
    
    local rankTable, i, A, r;
    
    rankTable := table([seq(i=0, i=0..n)]);
    
    for i to nops(S)^(n^2) do
        A := getMatrixAtIndex(i, n, S);
        r := LinearAlgebra:-Rank(A);
        rankTable[r] := rankTable[r] + 1;
    end do;
    
    return rankTable:
    
end proc: