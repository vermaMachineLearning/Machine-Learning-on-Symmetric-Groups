# Parameters
#    N::Int
#    - the problem size
# Return Values
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1} (Young's Orthogonal Representations)
#    - YOR[n][p][k] is Young's Orthogonal Representation for the Adjacent Transposition (K, K + 1) for the pth Partition of n
#    PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#    - output1 from partition_tree()
def yor(N):
    P, WI = partitions(N)
    PT = partition_tree(N, P, WI)
    YS = ys_symbols(N, P, PT)
    YOR = np.empty(shape=(N,), dtype=object)
    for n in range(1, N+1):
        YSn = YS[n-1]
        Pn = P[n-1]
        Pn_L = len(Pn)
        YORn = np.empty(shape=(Pn_L,), dtype=object)
        for p in range(1, Pn_L+1):
            YSnp = YSn[p-1]
            Pnp = Pn[p-1]
            R = len(Pnp)
            YORn[p-1] = yor_p(n, R, YSnp)
        
        YOR[n-1] = YORn
    
    return YOR, PT


# Parameters:
#    N::Int
#    - the size of P
#    R::Int
#    - the length of P
#    P::Array{Int, 1}
#    - P is a Partition of N
# Return Values:
#    YORp::Array{SparseMatrixCSC, 1}
#    - YORp[k] is the Young's Orthogonal Representation of P for the Adjacent Transposistion (k, k + 1)
def yor_p(N, R, YSymbols):
    YORp = np.empty(shape=(N-1,), dtype=object)
    L = len(YSymbols)
    DA, IA = ys_information(N, R, YSymbols)
    for k in range(1, N):
        YORp[k-1] = yor_pk(DA, IA, L, k)
    
    return YORp


#Parameters
#    DA::Array{Int, 2}  
#    - output1 from ys_information()
#    IA::Array{Int, 2}
#    - outpu2 from ys_information()
#    L::Int
#    - the number of Yamanouchi Symbols for the current Partition
#    K::Int
#    - represents the Adjacent Transposition (K, K + 1)
#Return Values
#    YORpk::SparseMatrixCSC 
#    - Young's Orthongal Representation of the Adjacent Transposition (K, K + 1) for the Partition p
def yor_pk(DA, IA, L, k):
    n = L
    for i in range(1, L+1):
        if IA[i-1,k-1] != 0:
            n += 1
        
    
    colptr = np.zeros(shape=(L+1,), dtype=np.int32)
    colptr[0] = 1
    rowval = np.zeros(shape=(n,), dtype=np.int32)
    nzval = np.zeros(shape=(n,), dtype=np.float32)
    cp = 1
    for i in range(1, L+1):
        D = 1/DA[i-1, k-1]
        I = IA[i-1,k-1]
        if I == 0:
            colptr[i] = colptr[i-1] + 1
            rowval[cp-1] = i
            nzval[cp-1] = D
            cp += 1
        else:
            colptr[i] = colptr[i-1] + 2
            if I > i:
                rowval[cp-1] = i
                nzval[cp-1] = D
                cp += 1
                rowval[cp-1] = I
                nzval[cp-1] = np.sqrt(1 - D * D)
                cp += 1
            else:
                rowval[cp-1] = I
                nzval[cp-1] = np.sqrt(1 - D * D)
                cp += 1
                rowval[cp-1] = i
                nzval[cp-1] = D
                cp += 1
            
    YORpk = scipy.sparse.csc_matrix((nzval,rowval-1,colptr-1),shape=(L, L))
    return YORpk


