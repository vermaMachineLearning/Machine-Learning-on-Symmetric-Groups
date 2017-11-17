# Parameters:
#    N::Int
#    - the problem size
#    K::Int
#    - the problem is homogenous at N-K
#    FFT::Array{Array{Float64, 2}, 1}
#    - a Fast Fourier Transform of size N
#    - should be the output of sn_fft_bl()
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor_bl()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor_bl()
#    ZFI::Array{Int, 1}
#    - output3 from yor_bl()
# Return Values:
#    SNF::Array{Float64, 1}
#    - the bandlimited function over Sn that corresponds to FFT
def sn_ifft_bl(N, K, FFT, YOR, PT, ZFI):
    
    SNF = np.zeros(shape=(factorial_ratio(N,K),), dtype=np.float32)
    compute_ifft_bl(N, K, SNF, FFT, YOR, PT, ZFI, Counter(1))
    return SNF

def factorial_ratio(N,K):
	ratio=1
	for i in range(N, N-K, -1):
		ratio=ratio*i
	return ratio


# Parameters:
#    N::Int
#    - the problem size
#    K::Int
#    - the problem is homogenous at N-K
#    SNF::Array{Float64, 1}
#    - the bandlimited function over Sn that is being calculated
#    FFT::Array{Array{Float64, 2}, 1}
#    - a Fast Fourier Transform of size N
#    - should be the output of sn_fft_bl()
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor_bl()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor_bl()
#    ZFI::Array{Int, 1}
#    - output3 from yor_bl()
#    C::Counter
#    - a wrapper for an Int that is used to coordinate the recursive defining of the elements of SNF
# Return Values:
#    None
#    - once the recursion is finished, SNF is ful
def compute_ifft_bl(N, K, SNF, FFT, YOR, PT, ZFI, C):
    
    if K == 0:
        V = FFT[-1][0, 0]
        V /= np.math.factorial(N)
        SNF[C.N-1] = V
        C.N += 1
    else:
        YORn = YOR[N-1]
        NPn = len(YORn)
        YORd = YOR[N-2]
        NPd = len(YORd)
        PTn = PT[N-1]
        ZFIn = ZFI[N-1]
        ZFId = ZFI[N-2]
        sFFT = np.empty(shape=(NPd,), dtype=object)
        
        for p in range(1, ZFId+1):
            sFFT[p-1] = np.zeros(shape=(1, 1), dtype=np.float32)
        
        for n in range(1, N+1):
            for p in range(ZFId + 1, NPd+1):
                Dim = YORd[p-1][0].shape[0]
                sFFT[p-1] = np.zeros(shape=(Dim, Dim), dtype=np.float32)
            
            for p in range(ZFIn + 1, NPn+1):
                update_sfft_bl(N, n, sFFT, FFT[p-1], YORn[p-1], YORd, PTn[p-1], ZFId)
            
            compute_ifft_bl(N-1, K-1, SNF, sFFT, YOR, PT, ZFI, C)
        
    


# Parameters:
#    N::Int
#    - the problem size
#    n::Int
#    - determines which left-sided coset we are defining the FFT for
#    sFFT::Array{Array{Float64, 2}, 1}
#    - the FFT of the nth left-sided coset that we are defining
#    FFTp::Array{Float64, 2}
#    - the pth component of the FFT of size N
#    YORnp::Array{SparseMatrixCSC, 1}, 1}
#    - Youngs Orthogonal Representations for the pth Partition of N as defined by yor_bl()
#    YORd::Array{Array{SparseMatrixCSC, 1}, 1}
#    - Youngs Orthogonal Reprsentations for the Partitions of N - 1 as defined by yor_bl()
#    PTnp::Array{Array{Array{Int, 1}, 1}, 1}
#    - The decomposition indices for the pth Partition of N as defined by yor_bl()
#    ZFId::Int
#    - if p <= ZFId, sFFT[p] is a zero-frequency component
# Return Values:
#    None
#    - sFFT is fully defined once all of the components of FFT have been used
def update_sfft_bl(N, n, sFFT, FFTp, YORnp, YORd, PTnp, ZFId):

    Dim = YORnp[0].shape[0]
    M = np.eye(Dim, dtype=np.float32)
    
    for ccn in range(n, N):
        M = np.matmul(M, YORnp[ccn-1].todense())
    
    M = np.transpose(M)
    M = np.matmul(M, FFTp)
    lb = 1
    
    for d in range(1, len(PTnp)+1):
        index = PTnp[d-1]
        if index > ZFId:
            sDim = YORd[index-1][0].shape[0]
            ub = lb + sDim - 1
            lb_pyidx=lb-1
            ub_pyidx=ub
            sFFT[index-1] += (Dim / (sDim * N))*M[lb_pyidx:ub_pyidx, lb_pyidx:ub_pyidx]
            lb += sDim
        else:
            lb += int(round(YORd[index-1][0][0, 0]))
			
    


# Parameters:
#    N::Int
#    - the problem size
#    K::Int
#    - the problem is homogenous at N-K
#    n::Int
#    - determines which left-sided coset we are defining the function for
#    FFT::Array{Array{Float64, 2}, 1}
#    - a Fast Fourier Transform of size N
#    - should be the output of sn_fft_bl()
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor_bl()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor_bl()
#    ZFI::Array{Int, 1}
#    - output3 from yor_bl()
# Return Values:
#    -SNF::Array{Float64, 1}
#    - the bandlimited function on the nth left-sided coset of Sn
def compute_sifft_bl(N, K, n, FFT, YOR,  PT, ZFI):
    
    pSNF = np.zeros(shape=(factorial_ratio(N,K),), dtype=np.float32)
    C = Counter(1)
    YORn = YOR[N-1]
    NPn = len(YORn)
    YORd = YOR[N-2]
    NPd = len(YORd)
    PTn = PT[N-1]
    ZFIn = ZFI[N-1]
    ZFId = ZFI[N-2]
    sFFT = np.empty(shape=(NPd,), dtype=object)
    
    for p in range(1, ZFId+1):
        sFFT[p-1] = np.zeros(shape=(1, 1), dtype=np.float32)
    
    for p in range(ZFId + 1, NPd+1):
        Dim = YORd[p-1][0].shape[0]
        sFFT[p-1] = np.zeros(shape=(Dim, Dim), dtype=np.float32)
    
    for p in range(ZFIn + 1, NPn+1):
        update_sfft_bl(N, n, sFFT, FFT[p-1], YORn[p-1], YORd, PTn[p-1], ZFId)
    
    compute_ifft_bl(N-1, K-1, pSNF, sFFT, YOR, PT, ZFI, C)
    return pSNF



