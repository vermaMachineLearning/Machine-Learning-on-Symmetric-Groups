# Parameters:
#    N::Int
#    - the problem size
#    SNF::Array{Foat64, 1}
#    - SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor()
# Return Values:
#    FFT::Array{Float64, 2}
#    - FFT is the Fast Fourier Transform of SNF
def sn_fft(N, SNF, YOR, PT):
    C = Counter(1)
    return compute_fft(N, SNF, YOR, PT, C) 

	
	# Parameters:
#    N::Int
#    - the size of the FFT being calculated is N
#    SNF::Array{Foat64, 1}
#    - SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1() from yor()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2() from yor()
#    C::Counter
#    - a wrapper for an Int that is used to coordinate the recursive access to the elements of SNF
# Return Values:
#    - FFT is the Fast Fourier Transfrom of SNF
def compute_fft(N, SNF, YOR, PT, C):
    
    if N == 1:
        sFFT = np.empty(shape=(1,), dtype=object)
        sFFTi = np.zeros(shape=(1, 1), dtype=np.float32)
        sFFTi[0,0] = SNF[C.N-1]
        sFFT[0] = sFFTi
        C.N += 1
        return sFFT
    
    sFFT = np.empty(shape=(N,), dtype=object)
    for n in range(1, N+1):
        sFFT[n-1] = compute_fft(N-1, SNF, YOR, PT, C)
    
    FFT = combine_sfft(N, YOR[N-1], PT[N-1], sFFT)
    return FFT



# Parameters:
#    N::Int
#    - the size of the FFT being calculated is N
#    YORn::Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - Young's Orthogonal Representations for the Partitions of N
#    PTn::Array{Array{Int, 1}, 1}
#    - the decomposition indices for the Partitions of N
#    sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#    - the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#    FFT::Array{Array{Float64, 2}, 1}
#    - FFT is the Fast Fourier Transform for this group
def combine_sfft(N, YORn, PTn, sFFT):
    
    NP = len(YORn)
    FFT = np.empty(shape=(NP,), dtype=object)
    
    for p in range(1, NP+1):
        YORnp = YORn[p-1]
        PTnp = PTn[p-1]
        FC = fc(N, YORnp, PTnp, sFFT)
        FFT[p-1] = FC
    return FFT


# Parameters:
#    N::Int
#    - the Fourier Component that is being calculated corresponds to a Partition, P, of N
#    YORnp::Array{SparseMatrixCSC, 1}
#    - Young's Orthogonal Representations for P
#    PTnp::Array{Int, 1}
#    - the indices for the Partitions that P decomposes into
#    sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#    - the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#    FC::Array{Float64, 2}
#    - FC is the Fourier Coefficient corresponding to the Partition P
def fc(N, YORnp, PTnp, sFFT):

    Dim = YORnp[0].shape[0]
    FC = dsm(Dim, sFFT[N-1], PTnp)
    CCM = np.eye(Dim)
    
    for n in range(N-1, 0, -1):
        CCM = np.matmul(YORnp[n-1].todense(), CCM)
        DSM = dsm(Dim, sFFT[n-1], PTnp)
        FC_n = np.matmul(CCM, DSM)
        FC += FC_n
        
    return FC

# Parameters:
#    Dim::Int
#    - the size of the DSM that is being caculated
#    sFFTn::Array{Array{Float64, 2}, 1}
#    - one of the elements of sFFT from fc()
#    PTnp::Array{Int, 1}
#    - same as in fc()
# Return Values:
#    DSM::Array{Float64, 2} (Direct Sum Matrix) 
#    - the direct sum of the coefficients of sFFTn that correspond to the Partitions indicated by PTnp 
def dsm(Dim, sFFTn, PTnp):
    
    DSM = np.zeros(shape=(Dim, Dim), dtype=np.float32)
    offset = 0
    
    for i in range(1, len(PTnp)+1):
        FC = sFFTn[PTnp[i-1]-1]
        sDim = FC.shape[0]
        for r in range(1, sDim+1):
            for c in range(1, sDim+1):
                v = FC[r-1, c-1]
                if v != 0.0:
                    DSM[offset+r-1, offset+c-1] = v
        offset += sDim
    return DSM


