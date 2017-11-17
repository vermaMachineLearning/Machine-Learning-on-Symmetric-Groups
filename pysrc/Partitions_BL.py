# Parameters:
#   N::Int
#   - the problem size
#   K::Int
#   - the problem is homogenous at N-K
# Return Values:
#   P::Array{Array{Array{Int, 1}, 1}, 1} (Partitions)
#   - P[n][p] contains the pth Partition of n that is required for the bandlimited functionality
#   - This is equivalent to the Partitions that correspond to non-zero frequencies in FFT and their immediate descendents
#   - length(P[n]) = 0 for n = 1:(N-K-1)
#   ZFI::Array{Int, 1} (Zero Frequency Information)
#   - ZFI[n] = k if, for p<=k, P[n][p] is a zero frequency partition
#   WI::Array{Int, 2} (Width Information)
#   - WI[n, w] contains the number of Partitions of n whose first element is less than or equal to w
#   - This only counts the Partitions that are in P and not all Partitions
def partitions_bl(N, K):

    P_K, WI_K = partitions(K)
    WI = np.zeros(shape=(N, N),dtype=np.int32)

    for n in range(N-K, N+1):
        WI[n-1,n-1] = 1

    for n in range(N-K, N):
        for w in range(N-K-1, n):
            WI[n-1,w-1] = WI_K[n-w-1, np.min((n-w, w))-1]


    for w in range(N-K, N):
        WI[N-1,w-1] = WI_K[N-w-1, np.min((N-w, w))-1]

    for n in range(N-K, N+1):
        for w in range(N-K, n+1):
            WI[n-1, w-1] += WI[n-1, w-2]

    P = np.empty(shape=(N,), dtype=object)
    ZFI = np.zeros(shape=(N,), dtype=np.int32)

    for n in range(1, N-K):
        P[n-1] = np.empty(shape=(0,), dtype=object)
        ZFI[n-1] = 0

    for n in range(N-K, N):
        Pn = np.empty(shape=(WI[n-1,n-1],), dtype=object)
        i = 1
        for spi in range(1, WI_K[n-N+K, np.min((n-N+K+1, N-K-1))-1]+1):
            Pn[i-1] = np.concatenate(([N-K-1], P_K[n-N+K][spi-1]),axis=0)
            i += 1

        ZFI[n-1] = i - 1
        for w in range(N-K, n):
            for spi in range(1, WI_K[n-w-1, np.min((n-w, w))-1]+1):
                Pn[i-1] = np.concatenate(([w],P_K[n-w-1][spi-1]))
                i += 1

        Pn[i-1] = np.array([n])
        P[n-1] = Pn

    Pn = np.empty(shape=(WI[N-1,N-1],), dtype=object)
    i = 1

    for w in range(N-K, N):
        for spi in range(1, WI_K[N-w-1, np.min((N-w, w))-1]+1):
            Pn[i-1] = np.concatenate(([w], P_K[N-w-1][spi-1]))
            i += 1

    Pn[i-1] = np.array([N])
    P[N-1] = Pn
    ZFI[N-1] = 0
    return P, ZFI, WI


