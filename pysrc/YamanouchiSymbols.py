# Let X be a Yamanouchi Symbol corresponding to a Partition of N whose length is R
#    X::Array{Int, 1}
#    length(X) = N
#    maximum(X) = R
#    X[i] = r means that i appears in row r of the Standard Tableau represented by X

# Parameters:
#    N::Int
#    - the problem size
#    P::Array{Array{Array{Int, 1}, 1}, 1}
#    - output1 from partitions()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output1 from partition_tree()
# Return Values
#    YS::Array{Array{Array{Array{Int, 1}, 1}, 1}, 1} (Yamanouchi Symbols)
#    - YS[n][p][i] is the ith Yamanouchi Symbol of the pth Partition of n
def ys_symbols(N, P, PT):
    
    YS = np.empty(shape=(N,), dtype=object)
    YSn = np.empty(shape=(1,), dtype=object)
    YSnp = np.empty(shape=(1,), dtype=object)
    YSnp[0] = np.array([1])
    YSn[0] = YSnp
    YS[0] = YSn
    for n in range(2, N+1):
        Pn = P[n-1]
        PTn = PT[n-1]
        Pn_L = len(Pn)
        YSn = np.empty(shape=(Pn_L,), dtype=object)
        for p in range(1, Pn_L+1):
            Pnp = Pn[p-1]
            Pnp_L = len(Pnp)
            PTnp = PTn[p-1]
            PTnp_L = len(PTnp)
            c = 1
            row = 1
            YSnp = np.empty(shape=(degree_v2(n, Pnp),), dtype=object)
            for d in range(1, PTnp_L+1):
                while row < Pnp_L and Pnp[row-1] <= Pnp[row]:
                    row += 1
                
                YSd = YS[n-2][PTnp[d-1]-1]
                YSd_L = len(YSd)
                for i in range(1, YSd_L+1):
                    YSnpc = np.zeros(shape=(n,), dtype=np.int32)
                    YSdi = YSd[i-1]
                    for j in range(1, n):
                        YSnpc[j-1] = YSdi[j-1]
                    
                    YSnpc[n-1] = row
                    YSnp[c-1] = YSnpc 
                    c += 1
                
                row += 1
            
            YSn[p-1] = YSnp
        
        YS[n-1] = YSn
    
    return YS


# Parameters:
#    N::Int
#    - the size of the Partition that the elements of YSymbols correspond to
#    R::Int  
#    - the length of the Partition that the elements of YSymbols correspond to
#    YSymbols::Array{Array{Int, 1}, 1}
#    - the array of Yamanouchi Symbols corresponding to a particular Partition
# Return Values
#    DA::Array{Int, 2} (Distance Array)
#    - output1 from ys_distance()
#    IA::Array{Int, 2} (Index Array)
#    - output1 from ys_indices()
def ys_information(N, R, YSymbols):
    CA = ys_content(N, R, YSymbols)
    DA = ys_distance(CA)
    IA = ys_indices(N, YSymbols, DA)
    return DA, IA


#Parameters:
#    N::Int
#    - the size of the Partition that the elements of YSymbols correspond to
#    R::Int  
#    - the length of the Partition that the elements of YSymbols correspond to
#    YSymbols::Array{Array{Int, 1}, 1}
#    - the array of Yamanouchi Symbols corresponding to a particular Partition
#Return Value
#    CA::Array{Int, 2} (Content Array)
#    - CA[i, n] is the Content of n for the ith Standard Tableau
def ys_content(N, R, YSymbols):
    
    CA = np.zeros(shape=(len(YSymbols),N), dtype=np.int32)
    
    for i in range(1, len(YSymbols)+1):
        YSymbol = YSymbols[i-1]
        PF = np.zeros(shape=(R,), dtype=np.int32)
        for n in range(1, N+1):
            row = YSymbol[n-1]
            PF[row-1] += 1
            CA[i-1,n-1] = PF[row-1] - row
    return CA


# Parameters
#    CA::Array{Int, 2}
#    - output1 from ys_content()
# Return Values
#    DA::Array{Int, 2} (Distance Array)
#    - DA[i, k] is the Distance of the ith Standard Tableau for the Adjacent Transposition (k, k + 1)
def ys_distance(CA):
    
    L, N = CA.shape
    DA = np.zeros(shape=(L, N-1), dtype=np.int32)
    for i in range(1, L+1):
        for n in range(1, N):
            DA[i-1, n-1] = CA[i-1, n] - CA[i-1, n-1]

    return DA


#Parameters
#    N::Int
#    - the size of the Partition that the elements of YSymbols correspond to
#    YSymbols::Array{Array{Int, 1}, 1}
#    - the array of Yamanouchi Symbols corresponding to a particular Partition
#    DA::Array{Int, 2}
#    - output1 from ys_distance()
#Return Values
#    IA::Array{Int, 2} (Index Array)
#    - IA[i, k] contains the index of the Standard Tableau that is the result of applying the Adjacent Transposition (k, k + 1) to the Standard Tableau at index i in YSymbols
#    - IA[i, k] = 0 when applying the Adjacent Transposition (K, K + 1) to the Standard Tableau at index i doesn't result in a Standard Tableau 
#Notes 
#    - there are block patterns that appear in this that could be used to speed up this search
#    - there may be ways to avoid searching all together
def ys_indices(N, YSymbols, DA):
    L = len(YSymbols)
    IA = np.zeros(shape=(L, N-1), dtype=np.int32)
    PF = np.zeros(shape=(L, N-1), dtype=np.int32)
    for i in range(1, L+1):
        YSymbol = YSymbols[i-1]
        k = 1
        ti = i + 1
        while k < N:
            D = DA[i-1, k-1]
            if D != 1 and D != -1 and PF[i-1, k-1] == 0:
                while ti <= L:
                    if istranspose(YSymbol, YSymbols[ti-1], k) == 1:
                        IA[i-1, k-1] = ti
                        IA[ti-1, k-1] = i
                        PF[i-1, k-1] = 1
                        PF[ti-1, k-1] = 1
                        break
                    ti += 1
            k += 1
    return IA


# Parameters
#    YSymbol1::Array{Int, 1}
#    - the starting Yamanouchi Symbol
#    YSmobl2::Array{Int, 2}
#    - the potential resulting Yamanouchi Symbol
#    K::Int
#    - represents the Adjacent Transposition (K, K + 1)
# Return Values
#    - 1, if YSymbol2 is the result of applying the Adjacent Transposition (K, K + 1) to YSymbol1
#    - 0, otherwise
def istranspose(YSymbol1, YSymbol2, K):
    
    for i in range(1, K):
        if YSymbol1[i-1] != YSymbol2[i-1]:
            return 0
        
    if YSymbol1[K-1] != YSymbol2[K] or YSymbol1[K] != YSymbol2[K-1]:
        return 0
    
    for i in range(K+2, len(YSymbol1)+1):
        if YSymbol1[i-1] != YSymbol2[i-1]:
            return 0
        
    return 1


