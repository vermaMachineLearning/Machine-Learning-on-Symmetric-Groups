# Parameters:
#   N::Int
#   - the problem size
#   PA::Array{Array{Int, 1}, 1}
#   - PA[i] is a Permutation of N
#   VA::Array{Float64, 1}
#   - VA[i] is the Value associated with PA[i]
# Return Values:
#   SNF::Array{Float64, 1}
#   - SNF[i] is the value associated with the permutation that permutation_index() maps to i
#   - this is the format for the SNF parameter of sn_fft()
# Notes:
#   - any permutation of N not represented in PA will be assigned a value of zero
def snf(N, PA, VA):
    SNF = np.zeros(np.math.factorial(N), dtype=np.float32)
    for i in range(1, len(PA)+1):
        Permutation = PA[i-1]
        Index = permutation_index(Permutation)
        SNF[Index-1] = VA[i-1]
    return SNF


# Parameters:
#   N::Int
#   - the problem size
#   K::Int
#   - the problem is homogenous at N-K
#   PA::Array{Array{Int, 1}, 1}
#   - PA[i] is a Permutation of N
#   VA::Array{Float64, 1}
#   - VA[i] is the Value associated with PA[i]
# Return Values:
#   SNF::Array{Float64, 1}
#   - SNF[i] is the value associated with all of the permutations that permutation_index() maps to any value in the range ((i - 1) * factorial(N - K) + 1):(i * factorial(N - K))
#   - this is the format for the SNF parameter of sn_fft_bl()
# Notes:
#   - any homogenous coset that doesn't have a representative permutation in PA will be assigned a value of zero
def snf_bl(N, K, PA, VA):
    #N_F = np.math.factorial(N)
    BS = np.math.factorial(N-K)
    snf_len=1
    SNF = np.zeros(np.round(factorial_ratio(N,K)).astype(np.int32)) #Python: Fix overflow problem!
    for i in range(1, len(PA)+1):
        Permutation = PA[i-1]
        Index = permutation_index(P)
        Index = np.ceil(Index / BS)
        SNF[Index-1] = VA[i-1]
    return SNF


# Parameters:
#   N::Int
#   - the problem size
#   PA::Array{Array{Int, 1}, 1}
#   - PA[i] is a Permutation of N
#   VA::Array{Float64, 1}
#   - VA[i] is the Value associated with PA[i]
# Return Values:
#   SNF::Array{Float64, 1}
#   - SNF[i] is the value associated with the permutation that permutation_index() maps to NZL[i]
#   - this is the format for the SNF parameter of sn_fft_sp()
#   NZL::Array{Int, 1}
#   - NZL must in increasing order
#   - this is the format for the NZL parameter of sn_fft_sp()
# Notes:
#   - the values in VA should be non-zero
def snf_sp(N, PA, VA):
    L = len(PA)
    SNF = np.zeros(L, dtype=np.float32)
    NZL = np.zeros(L, dtype=np.int32)
    SNF[0] = VA[0]
    NZL[0] = permutation_index(PA[0])
    for i in range(2, L+1):
        val = VA[i-1]
        index = permutation_index(PA[i-1])
        j = 0
        while j <= (i - 1):
            if index < NZL[j-1]:
                break
            j += 1
        for k in range(i, j,-1):
            SNF[k-1] = SNF[k-2]
            NZL[k-1] = NZL[k-2]

        SNF[j-1] = val
        NZL[j-1] = index

    return SNF, NZL

