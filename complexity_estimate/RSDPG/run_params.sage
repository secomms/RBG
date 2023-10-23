#########################
# check R-SDP(G) params #
#########################

def rts(t):
    s = str(round(t,1))
    return (s + "0" * (2 + s.find(".") -len(s)))


def stern(q, z, n, k, m):
    """
    inputs:
    p:      field size
    z:      size of restricted set
    n:      code length
    k:      code dimension
    m:      order of group G

    output:
    time complexity and parameters of optimized stern solver
    """


    C_opt = 10**9
    Cstern = 0

    Z = log(z,2)
    nsol = N(z**m * q**(k-n))+1

    for l in range(n-k):
            # number of enumerated error vectors
            La = z**ceil((k+l)/2)
            Lb = z**floor((k+l)/2)

            # cost of collision search
            Ccoll = La * Lb * q**(-l) 

            # overall cost
            Cstern = N(log((La + Lb + Ccoll)/nsol ,2)) #
            if Cstern < C_opt :
                C_opt = Cstern
                P_opt = {'ell': l}
    return C_opt, P_opt


def get_num_subcodes(z, n, m, w, d): 
    """
    get the expected number of subcodes with given parameters

    inputs:
    z: field size (here z!)
    n: length of C
    m: dimension of code (here m)
    w: support size of subcode
    d: dimension of subcode

    output:
    expected number of subcodes with given parameters
    """
    return binomial(n,w) *(z**d-_sage_const_1 )**(w-d) * q_binomial(m,d,z) / q_binomial(n,d,z)

# Table 2
qs = [1019, 347, 719, 971, 643, 269, 227, 107, 83, 223, 103, 79, 59, 47, 23, 53]
zs = [ 509, 173, 359,  97, 107,  67, 113,  53, 41,  37,  17, 13, 29, 23, 11, 13]
ns = [  40,  41,  49,  44,  60,  52,  43,  53, 73,  56,  76, 82, 63, 69, 93, 82]
ks = [  16,  20,  17,  26,  25,  27,  22,  26, 28,  33,  44, 49, 31, 34, 46, 47]
ms = [  18,  23,  20,  26,  26,  29, 24,   31, 35,  34,  48, 54, 38, 42, 61, 54]

# table 3
qs += [ 53, 103, 223, 1019] 
zs += [ 13,  17,  37,  509]
ns += [ 82,  76,  56,   40] 
ks += [ 47,  44,  33,   16] 
ms += [ 54,  48,  34,   18]  

seclvl = 128

for i in range(len(qs)):
    q = qs[i]
    n = ns[i]
    k = ks[i]
    R = k/n
    m = ms[i]
    z = zs[i]
    Z = log(z,2)
    Q = log(q,2)
    print("for q=",q,", z=",z," n=",n,"k=",k,"m=",m,)

    relevant_submatrix = 0
    # submatrices with width <0 m
    for j in range(1,m+1):
        for rho in range(1, j+1):
            if rho <= seclvl / Z:
                numcodes = get_num_subcodes(z,n,n-m,j,j-rho)
                if numcodes >= 1.0:
                    relevant_submatrix += 1


    # submatrices with width > m
    for j in range(m+1, n+1):
        for rho in range(1, m+1):
            if rho <= seclvl / Z:
                numcodes = get_num_subcodes(z,n,m,n-j,m-rho)
                if numcodes >= 1.0:
                    relevant_submatrix += 1

    C_stern, P_stern = stern(q, z, n, k, m)

    print(f'{relevant_submatrix} relevant submatrices per group; stern costs {rts(C_stern)} bit')
    print()
    
    