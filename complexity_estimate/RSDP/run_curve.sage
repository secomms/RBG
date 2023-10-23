load('evenz/BJMM2_even.sage')
load('oddz/BJMM2_odd.sage')
load('addStruct.sage')

shift   = True
onlyMax = True
avg     = False
if shift:   print("%shifted")
if onlyMax: print("%only best of difference set")

if avg:     print("%using avg values")
else:       print("%using max values")

num = 5
delta = 0.05
max_trials = 150
prec = 25

q = 157
z = 12
#z = 13
R = 0.45


Q = log(q,2)
Z = log(z,2)

print("%for (q=",q,", z=",z," R=",R,")")


E = get_E(q,z)

if shift:
    shift = E[0]
    E = shift_E(E, shift)

D = get_Diff(E, onlyMax)
zd = len(D)

nd = get_nd(E, True, avg)
nb = get_nb(E, D, True, avg)
print("% zd=",round_to_str(zd),"nd=",round_to_str(nd),"nb=",round_to_str(nb))



for w in range(1,prec):
    W = w/prec
    if hbin(W) + W*Z < (1-R)*Q and (not shift or W > 0.8):
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
                print("%    ", end = '' )
                for t in params:
                    i+=1
                    
                    print(t, round_to_str(params[t]),"    ", end = '' )
                    if i == 3:
                        i = 0
                        print()
                        print("%    ", end = '' )
                print()
                print("(",round_to_str(W),",",round_to_str(C),") %", round_to_str(C/Q))
                print("%---------------")   
