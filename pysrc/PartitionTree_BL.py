# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	P::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partitions_bl()
#	ZFI::Array{Int, 1}
#	- output2 from partitions_bl()
#	WI::Array{Int, 2}
#	- output3 from partitions_bl()
# Return Values:
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- for each value, i, in PT[n][j], P[n][j] decomposes into P[n-1][i]
#	- length(PT[n]) = 0 for n = 1:(N - K)
#	- length(PT[n][p]) = 0 for p <= ZFI[n]
def partition_tree_bl(N, K, P, ZFI, WI):
    
    PT = np.empty(shape=(N,), dtype=object)
    
    for n in range(1, N-K+1):
        PT[n-1] = np.empty(shape=(0,), dtype=object)
    
    for n in range(N-K+1, N+1):
        Pd = P[n-2]
        Pn = P[n-1]
        PTn = np.empty(shape=(len(Pn),), dtype=object)
        
        for p in range(1, ZFI[n-1]+1):
            PTn[p-1] = np.empty(shape=(0,), dtype=np.int32)
                
        for p in range(ZFI[n-1]+1, len(Pn)+1):
            
            PDA = pda(Pn[p-1])
            
            if (len(PDA) == 1 or PDA[0][0] == PDA[1][0]):
                lb = 1
                
                if (PDA[0][0] != 1):
                    lb = WI[n-2, PDA[0][0]-2]+1
                
                PTn[p-1] = dia(PDA, Pd, lb)
            else:
                lb1 = 1
                
                if (PDA[0][0] != 1):
                    lb1 = WI[n-2, PDA[0][0]- 2]+1
                
                lb2 = WI[n-2, PDA[1][0]-2]
                PTn[p-1] = dia_v2(PDA, Pd, lb1, lb2)
            
        
        PT[n-1] = PTn
    
    return PT


