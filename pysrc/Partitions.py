
def partitions(N):

    WI = np.zeros(shape=(N, N),dtype=np.int64)

    for n in range(1,N+1):
        WI[n-1, 0] = 1
    
    for w in range(2, N+1):
        WI[w-1, w-1] = 1
        for n in range(w+1, N+1):
            num = 0
            for i in range(1, np.min((n-w, w))+1):
                num += WI[n-w-1, i-1]
            WI[n-1, w-1] = num

    for n in range(2, N+1):
        for w in range(2, n+1):
            WI[n-1, w-1] += WI[n-1, w-2]

    P = np.empty(shape=(N,), dtype=object)
    Pn = np.empty(shape=(1,), dtype=object)
    Pn[0] = np.array([1])
    P[0] = Pn
    
    for n in range(2, N+1):
        Pn = np.empty(shape=(WI[n-1, n-1],), dtype=object)
        Pn[0] = np.ones(shape=(n,), dtype=np.int64)
        i = 2
        for w in range(2, n):
            for spi in range(1, WI[n-w-1, np.min((n-w, w))-1]+1):
                Pn[i-1] = np.concatenate(([w], P[n-w-1][spi-1]),axis=0)
                i += 1
        Pn[i-1] = np.array([n])
        P[n-1] = Pn

    return P, WI