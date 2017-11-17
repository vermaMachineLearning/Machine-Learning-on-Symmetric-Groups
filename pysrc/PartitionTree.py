# Parameters:
#	N::Int
#	- the problem size
#	P::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partitions()
#	WI::Array{Int, 2}
#	- output2 from partitions()
# Return Values:
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- for each value, i, in PT[n][p], P[n][p] decomposes into P[n-1][i]
#	- length(PT[1]) = 0
def partition_tree(N, P, WI):
    
    PT = np.empty(shape=(N,), dtype=object)
    PT[0] = np.empty(shape=(0,), dtype=object)
    
    for n in range(2, N+1):
        Pd = P[n-2]
        Pn = P[n-1]
        PTn = np.empty(shape=(len(Pn),), dtype=object)
        
        for p in range(1, len(Pn)+1):
            PDA = pda(Pn[p-1])
            
            if (len(PDA) == 1 or  PDA[0][0] == PDA[1][0]):
                lb = 1
                if (PDA[0][0] != 1):
                    lb = WI[n-2, PDA[0][0]-2]+1
                
                PTn[p-1] = dia(PDA, Pd, lb)
            else:
                lb1 = 1
                if PDA[0][0] != 1:
                    lb1 = WI[n-2, PDA[0][0]-2]+1
                
                lb2 = WI[n-2, PDA[1][0]-2]
                PTn[p-1] = dia_v2(PDA, Pd, lb1, lb2)
            
        PT[n-1] = PTn
    return PT


# Parameters:
#	P::Array{Int, 1}
#	- P is a Partition
# Return Values:
#	PDA::Array{Array{Int, 1}, 1} (Partition Decomposition Array)
#	- PDA[i] is the the ith Partition that P decomposes into 
def pda(P):
    P_L = len(P)    
    num = 1
    for i in range(1, P_L):
        if P[i-1] > P[i]:
            num += 1
        
    PDA = np.empty(shape=(num,), dtype=object)
    c = 1
    for i in range(1, P_L):
        if P[i-1] > P[i]:
            D = copy.copy(P)
            D[i-1] -= 1
            PDA[c-1] = D
            c += 1
    
    if P[P_L-1] == 1:
        D = copy.copy(P[0:(P_L - 1)])
        PDA[c-1] = D
    else:
        D = copy.copy(P)
        D[P_L-1] -= 1
        PDA[c-1] = D

    return PDA


# Parameters:
#	PDA::Array{Array{Int, 1}, 1}
#	- output1 from pda()
#	Pd::Array{Array{Int, 1}, 1}
#	- Pd is array containing all the partitions of the same size as those in PDA
#	lb::Int
#	- the lower bound for the index where the first element of PDA will be found
# Return Values:
#	DIA::Array{Int, 1} (Decomposition Index Array)
#	- PDA[i] = Pd[DIA[i]]
def dia(PDA, Pd, lb):
    PDA_L = len(PDA)
    DIA = np.empty(shape=(PDA_L,), dtype=np.int32)
    c = 1
    p = lb
    while True:
        if np.array_equal(Pd[p-1],PDA[c-1]):
            DIA[c-1] = p
            c += 1
            if c > PDA_L:
                return DIA
        p += 1

# Parameters:
#	PDA::Array{Array{Int, 1}, 1}
#	- output1 from pda()
#	Pd::Array{Array{Int, 1}, 1}
#	- Pd is array containing all the partitions of the same size as those in PDA
#	lb1::Int
#	- the lower bound for the index where the first element of PDA will be found
#	lb2::Int
#	- the lower bound for the index where the second element of PDA will be found
# Return Values:
#	DIA::Array{Int, 1}
#	- PDA[i] = Pd[DIA[i]]
def dia_v2(PDA, Pd, lb1, lb2):
    PDA_L = len(PDA)
    DIA = np.empty(shape=(PDA_L,), dtype=np.int32)
    p = lb1
    while True:
        if np.array_equal(Pd[p-1],PDA[0]):
            DIA[0] = p
            break
        p += 1
    
    p = lb2
    c = 2
    while True:
        if np.array_equal(Pd[p-1],PDA[c-1]):
            DIA[c-1] = p
            c += 1
            if c > PDA_L:
                return DIA
        p += 1
		
