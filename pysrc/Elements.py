# Let X be a Permutation of N
#    X::Array{Int, 1}
#    length(X) == N
#    X[i] = j indicates that the item in position i is sent to position j

###
# Group Operations
###

# Parameters:
#    P1::Array{Int, 1}
#    - the first permutation
#    P2::Array{Int, 1}
#    - the second permutation
# Return Values:
#    Prod::Array{Int, 1}
#    - the permutation that is P1 * P2
# Notes:
#     - P1 and P2 must be permutations of the same size
def sn_multiply(P1, P2):
    N = len(P1)   
    Prod = np.zeros(shape=(N,), dtype=np.int32)
    for i in range(1, N+1):
        Prod[i-1] = P1[P2[i-1]-1]
    return Prod


# Parameters:
#    P::Array{Int, 1}
#    - a permutation
# Return Values:
#    Inv::Array{Int, 1}
#    - the permutation that is the inverse of P    
def sn_inverse(P):
    N = len(P)
    Inv = np.zeros(shape=(N,), dtype=np.int32)
    for i in range(1, N+1):
        Inv[P[i-1]-1] = i
    
    return Inv


###
# Permutation Constructors
###

# Parameters:
#    N::Int
#   - the size of the permutation
# Return Values:
#    P::Array{Int, 1}
#    - a random permutation of N
def sn_p(N):
    index = np.random.randint(low=1, high=np.math.factorial(N)+1)
    P = index_permutation(N, index)
    return P


# Parameters:
#    N::Int
#    - the size of the permutation
#    LB::Int
#    - the first position that is reassigned
#    UB::Int
#    - the last position that is reassigned
# Return Values:
#    CC::Array{Int, 1}
#    - the permutation of N that is the contiguous cycle [[LB, UB]]
#    - this is the permutation that sends LB to LB + 1, LB + 1 to LB + 2, ... ,  UB - 1 to UB, and UB to LB
# Notes
#    - 1 <= LB <= UB <= N
def sn_cc(N, LB, UB):
    CC = np.zeros(shape=(N,), dtype=np.int32)
    for i in range(1, LB):
        CC[i-1] = i
    
    for i in range(LB, UB):
        CC[i-1] = i + 1
    
    CC[UB-1] = LB
    for i in range(UB + 1, N+1):
        CC[i-1] = i
    
    return CC


# Parameters:
#    N::Int
#    - the size of the permutation
# Return Values:
#    CC::Array{Int, 1}
#    - a random contiguous cycle of N
def sn_cc(N):
    lb = np.random.randint(low=1, high=N+1)
    ub = np.random.randint(low=lb, high=N+1)
    CC = sn_cc(N, lb, ub)
    return CC


# Parameters:
#    N::Int
#    - the size of the permutation
#    K::Int
#    - the position that is being reassigned
# Return Values:
#    AT::Array{Int, 1}
#    - the permutation of N that is the adjacent transposition (K, K+1)
#    - this is the permutation that sends K to K + 1 and K + 1 to K
# Notes:
#    - 1 <= K < N
def sn_at(N, K):
    AT = np.zeros(shape=(N,), dtype=np.int32)
    for i in range(1, K):
        AT[i-1] = i
    
    AT[K-1] = K + 1
    AT[K] = K
    for i in range(K + 2, N+1):
        AT[i-1] = i
    
    return AT


# Parameters:
#    N::Int
#    - the size of the permutation
# Return Values:
#    AT::Array{Int, 1}
#    - a random adjacent transposition of N
def sn_at(N):
    k = np.random.randint(low=1, high=N)
    AT = sn_at(N, k)
    return AT


# Parameters:
#    N::Int
#    - the size of the permutation
#    I::Int
#    - the first postition that is being reassigned
#    J::Int
#    - the second position that is being reassigned
# Return Values:
#    Tr::Array{Int, 1}
#    - the permutation of N that is the transposition (I, J)
#    - this is the permutation that sends I to J and J to I
# Notes:
#    - 1 <= I <= N
#    - 1 <= J <= N
def sn_t(N, I, J):
    Tr = np.random.randint(low=1, high=N)
    for i in range(1, N+1):
        Tr[i-1] = i
    
    Tr[I-1] = J
    Tr[J-1] = I
    return Tr


# Parameters:
#    N::Int
#    - the size of the permutation
# Return Values:
#    Tr::Array{Int, 1}
#    - a random transposition of N
def sn_t(N):
    i = np.random.randint(low=1, high=N+1)
    j = np.random.randint(low=1, high=N+1)
    Tr = sn_t(N, i, j)
    return Tr


###
# Factorizations (and related operations) on the Left Coset Tree
###

# Parameters:
#    P::Array{Int, 1}
#    - a permutation
# Return Values:
#    CCF::Array{Int, 1}
#    - the Contiguous Cycle Factoriztion of P
#    - P = product for i = 1:(N - 1) of sn_cc(N, CCF[i], N + 1 - i)
def permutation_ccf(P):
    N = len(P)
    CCF = np.zeros(shape=(N-1,), dtype=np.int32)
    i = 1
    for j in range(N, 1, -1):
        CCF[i-1] = P[j-1]
        cc = sn_cc(N, P[j-1], j)
        cc_inv = sn_inverse(cc)
        P = sn_multiply(cc_inv, P)
        i += 1
    
    return CCF


# Parameters:
#    CCF::Array{Int, 1}
#    - a contiguous cycle factorization of some permutation
# Return Values:
#    Index::Int
#    - the unique index that the permutation corresponding to CCF maps to
def ccf_index(CCF):
    N = len(CCF) + 1
    Index = 1
    for i in range(1, N):
        N -= 1          
        if CCF[i-1] != 1:
            Index += (CCF[i-1] - 1) * np.math.factorial(N)
        
    
    return Index


# Parameters:
#    P::Array{Int, 1}
#    - a permutation
# Return Values:
#    Index::Int
#    - the unique index that P maps to
def permutation_index(P):
    ccf = permutation_ccf(P)
    index = ccf_index(ccf)
    return index


# Parameters:
#    N::Int
#    - the size of the permutation that maps to Index
#    Index::Int
#    - the index of some permutation of N
# Return Values:
#    CCF::Array{Int, 1}
#    - the contiguous cycle factorization that corresponds to the permutation that maps to Index
def index_ccf(N, Index):
    CCF = np.zeros(shape=(N-1,), dtype=np.int32)
    Index -= 1
    for i in range(1, N):
        q = np.floor(Index / np.math.factorial(N - i))
        Index -= q * np.math.factorial(N - i)
        CCF[i-1] = q + 1
    
    return CCF


# Parameters:
#    CCF::Array{Int, 1}
#    - a contiguous cycle factorization of some permutation
# Return Values:
#    P::Array{Int, 1}
#    - the permutation that corresponds to CCF
def ccf_permutation(CCF):
    N = len(CCF) + 1
    P = np.zeros(shape=(N,), dtype=np.int32)
    for i in range(1, N+1):
        P[i-1] = i
    
    for i in range(1, N):
        cc = sn_cc(N, CCF[i-1], N+1-i)
        P = sn_multiply(P, cc)
    
    return P


# Parameters:
#    N::Int
#    - the size of the permutation that maps to Index
#    Index::Int
#    - the index of some permutation of N
# Return Values:
#    P::Array{Int, 1}
#    - the permutation of N that maps to Index
def index_permutation(N, Index):
    ccf = index_ccf(N, Index)
    permutation = ccf_permutation(ccf)
    return permutation


# Parameters:
#    P::Array{Int, 1}
#    - a permutation
# Return Values
#    ATF::Array{Int, 1}
#    - the adjacent transposition factorization of P
#    - P = product for i = 1:length(ATF) of sn_at(N, ATF[i])
def permutation_atf(P):
    N = len(P)
    CCF = permutation_ccf(P)
    dATF = np.empty(shape=(len(CCF),), dtype=object)
    L = 0
    for i in range(1, len(CCF)+1):
        P = CCF[i-1]
        D = N-P
        P -= 1
        L += D
        N -= 1
        dATFs = np.zeros(shape=(D,), dtype=np.int32)
        for j in range(1, D+1):
            dATFs[j-1] = P + j
        
        dATF[i-1] = dATFs
    
    ATF = np.zeros(shape=(L,), dtype=np.int32)
    index = 1
    for i in range(1, len(CCF)+1):
        for j in range(1, len(dATF[i-1])+1):
            ATF[index-1] = dATF[i-1][j-1]
            index += 1
        
    
    return ATF


# Parameters:
#    P::Array{Int, 1}
#    - a permutation
#    YORnp::Array{SparseMatrixCSC, 1}
#    - YORnp[i] is Young's Orthogonal Representation for the adjacent transposition (i, i + 1) corresponding to the pth partition of n
# Return Values:
#     RM::Array{Float64, 2}
#    - Young's Orthogonal Representation of P corresponding to the pth partition of n
def yor_permutation(P, YORnp):
    ATF = permutation_atf(P)
    if len(ATF) == 0:
        Dim = YORnp[0].shape[0]
        RM = np.eye(Dim, dtype=np.float32)
        return RM
    else:
        RM = copy.copy(YORnp[ATF[0]-1].tondense())
        for i in range(2, len(ATF)+1):
            RM = RM * YORnp[ATF[i-1]-1].todense()
        
        return RM
    

