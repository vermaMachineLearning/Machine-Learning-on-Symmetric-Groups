# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
# Return Values
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1} (Young's Orthogonal Representations)
#	- YOR[n][p][k] is Young's Orthogonal Representation for the Adjacent Transposition (K, K + 1) for the pth Partition of n that is needed for the bandlimited functionality
#	- length(YOR[n]) = 0 for n = 1:(N - K - 1)
#	- if p < ZFI[n], length(YOR[n][p] = 1) and YOR[n][p][1,1] contains the dimension of the full Young's Orthogonal Representation
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 of partition_tree_bl()
#	ZFI::Array{Int, 1}
#	- output2 for paritions_bl()
def yor_bl(N, K):
	
	P, ZFI, WI = partitions_bl(N, K)
	PT = partition_tree_bl(N, K, P, ZFI, WI)
	YS = ys_symbols_bl(N, K, P, ZFI, PT)
	YOR = np.empty(shape=(N,), dtype=object)
	
	for n in range(1, N-K):
		YOR[n-1] = np.empty(shape=(0,), dtype=object)
	
	for n in range(N-K, N+1):
		YSn = YS[n-1]
		Pn = P[n-1]
		Pn_L = len(Pn)
		YORn = np.empty(shape=(Pn_L,), dtype=object)
		
		for p in range(1, ZFI[n-1]+1):
			YORnp_sparse = scipy.sparse.csc_matrix((1,1))
			YORnp_sparse[0,0] = len(YSn[p-1])
			YORnp = np.empty(shape=(1,), dtype=object)
			YORnp[0] = YORnp_sparse
			YORn[p-1] = YORnp
				
		for p in range(ZFI[n-1]+1, Pn_L+1):
			YSnp = YSn[p-1]
			Pnp = Pn[p-1]
			R = len(Pnp)
			YORn[p-1] = yor_p(n, R, YSnp)
		
		YOR[n-1] = YORn
	
	return YOR, PT, ZFI


