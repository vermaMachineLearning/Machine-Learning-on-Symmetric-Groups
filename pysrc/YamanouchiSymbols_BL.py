# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	P::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partitions_bl()
#	ZFI::Array{Int, 1}
#	- output2 from partitions_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partition_tree_bl()
# Return Values:
#	YS::Array{Array{Array{Array{Int, 1}, 1}, 1}, 1} (Yamanouchi Symbols)
#	- YS[n][p][i] is the ith Yamanouchi Symbol of the pth Partition of n is required for the bandlimited functionality
#	- length(YS[n]) = 0 for n = 1:(N - K - 1)
def ys_symbols_bl(N, K,  P, ZFI, PT):
    
    YS = np.empty(shape=(N,), dtype=object)
    
    for n in range(1, N-K):
        YS[n-1] = np.empty(shape=(0,), dtype=object)
    
    YSn = np.empty(shape=(len(P[N-K-1]),), dtype=object)
    
    for p in range(1, len(P[N-K-1])+1):
        Pnp = P[N-K-1][p-1]
        YSn[p-1] = ys_partition(N-K, Pnp)
    
    YS[N-K-1] = YSn
    
    for n in range(N-K+1, N+1):
        Pn = P[n-1]
        PTn = PT[n-1]
        Pn_L = len(Pn)
        YSn = np.empty(shape=(Pn_L,), dtype=object)
        for p in range(1, ZFI[n-1]+1):
            Pnp = Pn[p-1]
            YSn[p-1] = ys_partition(n, Pnp)
        
        for p in range(ZFI[n-1] + 1, Pn_L+1):
            Pnp = Pn[p-1]
            Pnp_L = len(Pnp)
            PTnp = PTn[p-1]
            PTnp_L = len(PTnp)
            c = 1
            row = 1
            YSnp = np.empty(shape=(degree_v2(n, Pnp),), dtype=object)
            for d in range(1, PTnp_L+1):
                while (row < Pnp_L and Pnp[row-1] <= Pnp[row]):
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


# Parameters
#	N::Int
#	- N is the size of P
#	P::Array{Int, 1}
#	- P is a Partition of N
# Return Values
#	YS_P::Array{Array{Int, 1}, 1} (Yamanouchi Sybol for the Partition P)
#	- YS_P[i] is the ith Yamanouchi Symbol of the Partition P
def ys_partition(N, P):
    if len(P) == 1:
        YS_P = np.empty(shape=(1,), dtype=object)
        YS_P[0] = np.ones(shape=(N,), dtype=np.int32)
        return YS_P
        
    NST = degree_v2(N, P)    
    YSymbol = np.zeros(shape=(N,), dtype=object)
    PA = copy.copy(P)
    YS_P = np.empty(shape=(NST,), dtype=object)
    C = Counter(NST)
    ys_fillin(N, YSymbol, PA, YS_P, C)
    return YS_P	


# Parameters:
#	n::Int
#	- the next element to be placed into YSymbol
#	YSymbol::Array{Int, 1}
#	- the partially completed Yamanouchi Symbol that is being filled in
#	PA::Array{Int, 1}
#	- PA[i] is the number of Positions Available in row i of the Standard Tableau corresponding to YSymbol
#	YS_P::Array{Array{Int, 1}, 1}
#	- the array that stores the Yamanouchi symbols as they are completed
#	C::Counter
#	- A wrapper for an int that allows the recursive calls to know the next empty index in YS_P
# Return Values:
#	NONE
#	- once the recursion is finished, YS_P is full 
# Notes:
#	- This shouldn't be called on its own, use its wrapper function ys_partition()
def ys_fillin(n, YSymbol, PA, YS_P, C):
    PA_L = len(PA)    
    if PA[PA_L-1] > 0:
        YSymbol_new = copy.copy(YSymbol)
        YSymbol_new[n-1] = PA_L
        PA_new = copy.copy(PA)
        PA_new[PA_L-1] -= 1
        ys_fillin(n-1, YSymbol_new, PA_new, YS_P, C)
    
    for i in range(PA_L-1, 1, -1): 
        if PA[i-1] > PA[i]:
            YSymbol_new = copy.copy(YSymbol)
            YSymbol_new[n-1] = i
            PA_new = copy.copy(PA)
            PA_new[i-1] -= 1
            ys_fillin(n-1, YSymbol_new, PA_new, YS_P, C)
        
    
    if PA[0] > PA[1]:
        if PA[0] == 1:
            YSymbol_new = copy.copy(YSymbol)
            YSymbol_new[0] = 1
            YS_P[C.N-1] = YSymbol_new
            C.N -= 1
        else:
            YSymbol_new = copy.copy(YSymbol)
            YSymbol_new[n-1] = 1
            PA_new = copy.copy(PA)
            PA_new[0] -= 1
            ys_fillin(n-1, YSymbol_new, PA_new, YS_P, C)
		
	


