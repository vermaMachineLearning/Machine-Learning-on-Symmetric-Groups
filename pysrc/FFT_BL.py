# Parameters:
#    N::Int
#    - the problem size N
#    K::Int
#    - the problem is homogenous at N-K
#    SNF::Array{Foat64, 1}
#    - SNF[i] is the value assigned to the ith homogenous subgroup of size N-K
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor_bl()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor_bl()
#    ZFI::Array{Int, 1}
#    - output3 from yor_bl()
# Return Values:
#    FFT::Array{Float64, 2}
#    - FFT is the Fast Fourier Transform of SNF
def sn_fft_bl(N, K, SNF, YOR, PT, ZFI):
    return compute_fft_bl(N, K, SNF, YOR, PT, ZFI, Counter(1))


# Parameters:
#    N::Int
#    - the problem size N
#    K::Int
#    - the problem is homogenous at N-K
#    SNF::Array{Foat64, 1}
#    - SNF[i] is the value assigned to the ith homogenous subgroup of size N-K
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor_bl()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor_bl()
#    ZFI::Array{Int, 1}
#    - output3 from yor_bl()
#    C::Counter()
#    - a wrapper for an Int that is used to coordinate the recursive access to the elements of SNF
# Return Values:
#    FFT::Array{Float64, 2}
#    - FFT is the Fast Fourier Transform of SNF
def compute_fft_bl(N, K, SNF, YOR, PT, ZFI, C):
    sFFT = np.empty(shape=(N,), dtype=object)
    if K != 1:
        for n in range(1, N+1):
            sFFT[n-1] = compute_fft_bl(N-1, K-1, SNF, YOR, PT, ZFI, C)
        
    else:
        for n in range(1, N+1):
            V = SNF[C.N-1]
            sFFT[n-1] = compute_fft_iv(N-1, V, YOR[N-2])
            C.N += 1 
    
    
    FFT = combine_sfft_bl(N, YOR[N-1], PT[N-1], sFFT, ZFI[N-1], ZFI[N-2])
    return FFT

# Parameters:
#    N::Int
#    - the size of the homogenous (invariant) subgroup
#    V::Float64
#    - the value associated with this homogenous subgroup
#    YORn::Array{SparseMatrixCSC, 1}
#    - Young's Orthogonal Representations for the Partitions of N as defined by yor_bl()
def compute_fft_iv(N, V, YORn):
    
    NP = len(YORn)
    FFT = np.empty(shape=(NP,), dtype=object)
    
    for p in range(1, NP):
        FFT[p-1] = YORn[p-1][0].todense()
    
    NZ = np.zeros(shape=(1,1),dtype=np.float32)
    NZ[0,0] = np.math.factorial(N)*V
    FFT[NP-1] = NZ
    
    return FFT


# Parameters:
#    N::Int
#    - the size of the FFT being calculated is N
#    YORn::Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - Young's Orthogonal Representations for the Partitions of N as defined by yor_bl()
#    PTn::Array{Array{Int, 1}, 1}
#    - the decomposition indices for the Partitions of N as defined by partition_tree_bl()
#    sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#    - the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#    FFT::Array{Array{Float64, 2}, 1}
#    - FFT is the Fast Fourier Transform for this group
def combine_sfft_bl(N, YORn, PTn, sFFT, ZFIn, ZFId):
    
    NP = len(YORn)
    FFT = np.empty(shape=(NP,), dtype=object)
    
    for p in range(1, ZFIn+1):
        FFT[p-1] = YORn[p-1][0].todense()
    
    for p in range(ZFIn + 1, NP+1):
        YORnp = YORn[p-1]
        PTnp = PTn[p-1]
        FC = fc_bl(N, YORnp, PTnp, sFFT, ZFId)
        FFT[p-1] = FC
    
    return FFT


# Parameters:
#    N::Int
#    - the Fourier Component that is being calculated corresponds to a Partition, P, of N
#    YORnp::Array{SparseMatrixCSC, 1}
#    - Young's Orthogonal Representations for P as defined by yor_bl()
#    PTnp::Array{Int, 1}
#    - the indices for the Partitions that P decomposes into as defined by pt_bl()
#    sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#    - the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#    FC::Array{Float64, 2}
#    - FC is the Fourier Coefficient corresponding to the Partition P
def fc_bl(N, YORnp, PTnp, sFFT, ZFId):
    
    Dim = YORnp[0].shape[0]
    FC = dsm_bl(Dim, sFFT[N-1], PTnp, ZFId)
    CCM = np.eye(Dim)
    
    for n in range(N-1, 0, -1):
        CCM = np.matmul(YORnp[n-1].todense(), CCM)
        DSM = dsm_bl(Dim, sFFT[n-1], PTnp, ZFId)
        FC_n = np.matmul(CCM, DSM)
        FC += FC_n
    
    return FC

# Parameters:
#    Dim::Int
#    - the size of the DSM that is being caculated
#    sFFTn::Array{Array{Float64, 2}, 1}
#    - one of the elements of sFFT from fc_bl()
#    PTnp::Array{Int, 1}
#    - same as in fc_bl()
# Return Values:
#    DSM::Array{Float64, 2} (Direct Sum Matrix)
#    - the direct sum of the coefficients of sFFTn that correspond to the Partitions indicated by PTnp 
def dsm_bl(Dim, sFFTn, PTnp, ZFId):
    
    DSM = np.zeros(shape=(Dim, Dim), dtype=np.float32)
    offset = 0
    
    for i in range(1, len(PTnp)+1):
        index = PTnp[i-1]
        if index > ZFId:
            FC = sFFTn[index-1]
            sDim = FC.shape[0]
            for r in range(1, sDim+1):
                for c in range(1, sDim+1):
                    v = FC[r-1, c-1]
                    if v != 0.0:
                        DSM[offset+r-1, offset+c-1] = v
            offset += sDim
        else:
            FC = sFFTn[index-1]
            offset += int(round(FC[0,0]))
    return DSM

