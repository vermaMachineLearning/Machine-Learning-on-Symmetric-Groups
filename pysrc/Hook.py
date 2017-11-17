#POTENTIAL ISSUE - line 22 may throw an Inexact Error for large problems

# Parameters:
#	N::Int
#	- N is the size of P
#	P::Array{Int, 1}
#	- P is a Partition
#Return Value:
#	NST::Int
#	- the number of Standard Tableau of the Young Diagram of P
def degree(N, P):
    print("N = ",N)
    print("P = ",P)
    if N > 10:
        NST = 3628800.0
        for i in range(1, len(P)+1):
            for j in range(1, P[i-1]):
                NST = NST/hook_length(P, i, j)
            
        
        for n in range(11, N+1):
            NST = NST*n
        
        NST = int(round(NST))
        print("NST = ",NST)
        return NST
    else:
        NST = np.math.factorial(N)
        NST = NST/hook_product(P)
        NST = int(round(NST))
        print("NST = ",NST)
        return NST

def log_factorial(N):
    
    log_factorial=0
    for i in range(1,N+1):
        log_factorial=log_factorial+np.log10(i)
    
    return log_factorial

def degree_v2(N, P):
    print("N = ",N)
    print("P = ",P)
        
    log_NST = log_factorial(N)-log_hook_product(P)
    NST=np.power(10,log_NST)
    NST = int(round(NST))
    print("NST = ",NST)
    return NST
	


# Parameters:
#	P::Array{Int, 1}
#	- P is a Partition
# Return Values:
#	HP::Int
#	- the Hookproduct of the Young Diagram of P
def log_hook_product(P):
    log_HP = 0
    for i in range(1, len(P)+1):
        for j in range(1, P[i-1]+1):
            log_HP = log_HP+np.log10(hook_length(P, i, j))
    return log_HP
   
def hook_product(P):
    HP = 1
    for i in range(1, len(P)+1):
        for j in range(1, P[i-1]+1):
            HP *= hook_length(P, i, j)
    return HP


# Parameters:
#	P::Array{Int, 1}
#	- P is a Partition
#	i::Int 
#	- the row index
#	j::Int
#	- the column index
# Return Value
#	HL::Int
#	- the HookLength of the box in the ith row and jth colomn of the Young Diagram of P
def hook_length(P, i, j):
    
    R = P[i-1] - j
    B = 0
    i += 1
    while (i <= len(P) and j <= P[i-1]):
        i += 1
        B += 1
    
    HL = 1+R+B
    return HL


