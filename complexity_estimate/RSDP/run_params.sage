load('BJMM2_even.sage')
load('BJMM2_odd.sage')
debug   = False
shift   = False
onlyMax = True
avg     = True
if shift:   print("shifted")
if onlyMax: print("only best of difference set")
if avg:     print("using avg values")

num = 5
delta = 0.5
max_trials = 150
# table 1
qs = [677, 379, 197, 2017, 1021, 197]
zs = [26,   21,  14,   63,  30,   14]
ns = [84,  103, 103,   70,  79,  102]
ks = [42,   52,  51,   32,  40,   51]
ws = [73,   82,  91,   70,  79,  102]

# table 3
qs += [ 67, 197, 991, 991] 
zs += [ 11,  14,  33,  33]
ns += [147, 105,  77,  77] 
ks += [ 63,  53,  48,  38] 
ws += [147, 105,  77,  77]  


for i in range(len(qs)):
    q = qs[i]
    n = ns[i]
    k = ks[i]
    R = k/n
    w = ws[i]
    W = w/n
    z = zs[i]
    Q = log(q,2)
    print("for (q=",q,", z=",z," n=",n,"k=",k,"w=",w,")")
    print("R=",round_to_str(R),"W=",round_to_str(W))
    
    if debug:
        E = get_E(q,z)
        if shift:
            shift = E[0]
            E = shift_E(E, shift)

    
        D = get_Diff(E, onlyMax)
        zd = len(D)

        nd = get_nd(E, True, avg)
        nb = get_nb(E, D, True, avg)
        print("% zd=",round_to_str(zd),"nd=",round_to_str(nd),"nb=",round_to_str(nb))
    C_best = 100
    for _ in range(num):
        if z%2 == 1 or shift:
            C, params, success = optimize_BJMM2odd(q, z, shift, onlyMax, avg, R, W, max_trials, False)
        else:
            C, params, success = optimize_BJMM2even(q, z, shift, onlyMax, avg, R, W, max_trials, False)

        if C < C_best-delta and success:
            C_best = C
            p_best = params
            i = 0
            if debug:
                print("%    ", end = '' )
                for t in params:
                    i+=1
                    
                    print(t, round_to_str(params[t]),"    ", end = '' )
                    if i == 3:
                        i = 0
                        print()
                        print("%    ", end = '' )
                print()
            print(round_to_str(C_best),"->",round_to_str(C_best*n))
            print("%---------------")   
