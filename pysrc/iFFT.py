# Parameters:
#    N::Int
#    - the problem size
#    FFT::Array{Array{Float64, 2}, 1}
#    - a Fast Fourier Transform of size N
#    - should be the output of sn_fft() or sn_fft_sp()
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor()
# Return Values:
#    SNF::Array{Float64, 1}
#    - the function over Sn that corresponds to FFT
def sn_ifft(N, FFT, YOR, PT):

    SNF = np.zeros(shape=(np.math.factorial(N), ), dtype=np.float64)
    compute_ifft(N, SNF, FFT, YOR, PT, Counter(1))
    return SNF
    

# Parameters:
#    N::Int
#    - the problem size
#    SNF::Array{Float64, 1}
#    - the function over Sn that is being calculated
#    FFT::Array{Array{Float64, 2}, 1}
#    - a Fast Fourier Transform of size N
#    - should be the output of sn_fft() or sn_fft_sp()
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor()
#    C::Counter
#    - a wrapper for an Int that is used to coordinate the recursive defining of the elements of SNF
# Return Values:
#    None
#    - once the recursion is finished, SNF is full
def compute_ifft(N, SNF, FFT, YOR, PT, C):
    
    if N == 2:
        YORn = YOR[1]
        sFFT = 0.5 * (YORn[0][0].todense() * FFT[0] + YORn[1][0].todense() * FFT[1])
        SNF[C.N-1] = sFFT[0,0]
        C.N += 1
        sFFT = 0.5 * (FFT[0] + FFT[1])
        SNF[C.N-1] = sFFT[0,0]
        C.N += 1
    else:
        YORn = YOR[N-1]
        NPn = len(YORn)
        YORd = YOR[N-2]
        NPd = len(YORd)
        PTn = PT[N-1]
        sFFT = np.empty(shape=(NPd, ), dtype=object)
        for n in range(1, N+1):
            for p in range(1, NPd+1):
                Dim = YORd[p-1][0].shape[0]
                sFFT[p-1] = np.zeros(shape=(Dim, Dim), dtype=np.float32)
            
            for p in range(1, NPn+1):
                update_sfft(N, n, sFFT, FFT[p-1], YORn[p-1], YORd, PTn[p-1])
            
            compute_ifft(N-1, SNF, sFFT, YOR, PT, C)
        
    


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
#    - Youngs Orthogonal Representations for the pth Partition of N
#    YORd::Array{Array{SparseMatrixCSC, 1}, 1}
#    - Youngs Orthogonal Reprsentations for the Partitions of N - 1
#    PTnp::Array{Array{Array{Int, 1}, 1}, 1}
#    - The decomposition indices for the pth Partition of N
# Return Values:
#    None
#    - sFFT is fully defined once all of the components of FFT have been used
def update_sfft(N, n, sFFT, FFTp, YORnp, YORd, PTnp):
    
    Dim = YORnp[0].shape[0]
    M = np.eye(Dim,dtype=np.float32)
    
    for ccn in range(n, N):
        M = np.matmul(M, YORnp[ccn-1].todense())
    
    M = np.transpose(M)
    M = np.matmul(M, FFTp)
    lb = 1
    
    for d in range(1, len(PTnp)+1):
        index = PTnp[d-1]
        sDim = YORd[index-1][0].shape[0]
        ub = lb + sDim - 1
        lb_pyidx=lb-1
        ub_pyidx=ub        
        sFFT[index-1] += (Dim / (sDim * N)) * M[lb_pyidx:ub_pyidx, lb_pyidx:ub_pyidx]
        lb += sDim
    


# Parameters:
#    N::Int
#    - the problem size
#    n::Int
#    - determines which left-sided coset we are defining the function for
#    FFT::Array{Array{Float64, 2}, 1}
#    - a Fast Fourier Transform of size N
#    - should be the output of sn_fft() or sn_fft_sp()
#    YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#    - output1 from yor()
#    PT::Array{Array{Array{Int, 1}, 1}, 1}
#    - output2 from yor()
# Return Values:
#    pSNF::Array{Float64, 1}
#    - the function on the nth left-sided coset of Sn
def compute_sifft(N, n, FFT, YOR,  PT):
    
    pSNF = np.zeros(shape=(np.math.factorial(N-1),), dtype=np.float32)
    C = Counter(1)
    YORn = YOR[N-1]
    NPn = len(YORn)
    YORd = YOR[N-2]
    NPd = len(YORd)
    PTn = PT[N-1]
    sFFT = np.empty(shape=(NPd,), dtype=object)
    
    for p in range(1, NPd+1):
        Dim = YORd[p-1][0].shape[0]
        sFFT[p-1] = np.zeros(shape=(Dim, Dim), dtype=np.float32)
    
    for p in range(1, NPn+1):
        update_sfft(N, n, sFFT, FFT[p-1], YORn[p-1], YORd, PTn[p-1])
    
    compute_ifft(N-1, pSNF, sFFT, YOR, PT, C)
    return pSNF


