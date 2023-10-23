"""
Optimization target: restricted BJMM2even

Run:
>>> optimize_BJMM2even()

"""
load('addStruct.sage')
load('macros.sage')
import collections
from math import*
import scipy.optimize as opt
import matplotlib.pyplot as plt
import numpy as np


set_BJMM2even = collections.namedtuple('BJMM2even', 'L p0 p1 p2 m1 m2 e1 e2 d1 d2 b1 b2 l0 l1 l2 l3 niter r0 r1')
def BJMM2even(f) : return wrap(f,set_BJMM2even)
#=================================
# Variables to be optimized:
# L:     redundancy of the small instance
# pi:    weight on level i
# li:    list size on lvl i
# niter: number of iteration
# ri:    representations for list on level i

scale = 3000

def time_BJMM2even(x):
    x = set_BJMM2even(*x)
    return (max(x.l3, x.l2, 2*x.l2 + x.r1 -x.r0, x.l0) +x.niter)*scale #scale up if no successful optimizazion

def BJMM2even_list(R,L,P,M,Z,ZD):
    """
    intermediate list sizes for BJMM2even
    R: rate
    L: red small instance
    P: small instance
    """
    return (R+L) * trinom(P/(R+L),M/(R+L)) + Z*P + ZD*M
    
def BJMM2even_repr(R,L, P, E, M, B, D, Z, ND, NB):
    """
    intermediate list sizes for BJMM2even
    R: rate
    L: red small instance
    pi: weight after merge
    ei: number of overlaps
    """
    repr = M
    t = R+L-P-M
    if t > 0.0:                         repr += t*hbin(E/t) + Z*E
    if P > 0.0 and 2*D + 2*B <= P:      repr += P + P*trinom(2*D/P, 2*B/P) + 2*ND*D + 2*NB*B
    return repr

def optimize_BJMM2even(q, z, shift, onlyMax, avg, R, W, max_trials, verb=True):
    """
    Optimizes the BJMM2even algorithm
    q: field size
    z: EVEN restricted size
    R: code rate
    W: weight of error
    """
    assert not shift,"shifting does not work with BJMM-even"
    assert z%2==0,"not even!"
    E = get_E(q,z)

    D = get_Diff(E, onlyMax)
    zd = len(D)

    nd = get_nd(E, verb, avg)
    nb = get_nb(E, D, verb, avg)
    if verb: print("zd=",zd,"nd=",nd,"nb=",nb)

    Q = log(q, 2)
    Z = log(z, 2)
    if verb and hbin(W) + W*Z > (1-R)*Q: print("uniqueness not given! badness=",hbin(W) + W*Z-(1-R)*Q)

    if zd > 0:  ZD = log(zd,2)
    else:       ZD = -1000

    if nd > 0:  ND = log(nd,2)
    else:       ND = -1000

    if nb > 0:  NB = log(nb,2)
    else:       NB = -1000



    time = time_BJMM2even
    objective = time
    mycons = [
    # restrictions
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : min(R+x.L, W)-x.p0)}, # P <= R+L
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : max(1-W-R-x.L,0)+x.p0)}, # P >= W-(1-R-L) = W + R + L -1
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : R + x.L - x.p0 - x.e1 - 0.0 )},
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : R + x.L - x.p1 - x.e2 - x.m1)},
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : x.p0/2 - x.d1 - x.m1)},
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : x.p1/2 - x.d2 - x.m2)},
    # parameters
    { 'type' : 'eq', 'fun' :   BJMM2even(lambda x : x.p0/2 + x.e1 + x.d1  - x.p1)}, 
    { 'type' : 'eq', 'fun' :   BJMM2even(lambda x : x.p1/2 + x.e2 + x.d2  - x.p2)},
    { 'type' : 'eq', 'fun' :   BJMM2even(lambda x :            x.b1 - x.m1)}, 
    { 'type' : 'eq', 'fun' :   BJMM2even(lambda x : (x.m1)/2 + x.b2 - x.m2)}, 
    # success probability
    { 'type' : 'eq', 'fun' :  BJMM2even(lambda x : Niter(R,W,x.L,x.p0) - x.niter)},
    # number of representations 
    { 'type' : 'eq', 'fun' :   BJMM2even(lambda x : BJMM2even_repr(R,x.L, x.p1, x.e2, x.m1, x.b2, x.d2, Z, ND, NB) - x.r1)},
    { 'type' : 'eq', 'fun' :   BJMM2even(lambda x : BJMM2even_repr(R,x.L, x.p0, x.e1, 0.0 , x.b1, x.d1, Z, ND, NB) - x.r0)},
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : Q*x.L - x.r0)}, # r0 <= Q*L
    { 'type' : 'ineq', 'fun' : BJMM2even(lambda x : x.r0  - x.r1)},
    # sizes of the lists
    { 'type' : 'eq', 'fun' : BJMM2even(lambda x : BJMM2even_list(R/2,x.L/2,x.p2/2, x.m2/2, Z, ZD)   - x.l3 )}, # base list 
    { 'type' : 'eq', 'fun' : BJMM2even(lambda x : 2*x.l3                                    - x.r1  - x.l2 )}, # concat-merged  
    { 'type' : 'eq', 'fun' : BJMM2even(lambda x : BJMM2even_list(R, x.L, x.p1, x.m1, Z, ZD) - x.r0  - x.l1 )}, # repr-merged
    { 'type' : 'eq', 'fun' : BJMM2even(lambda x : 2*x.l1                             + x.r0 - Q*x.L - x.l0 )}, # repr-merged 
    ]
    #         'L p0 p1 p2 m1 m2 e1 e2 d1 d2 b1 b2 l0 l1 l2 l3 niter r0 r1
    bounds = [(0, 1-R)]+[(max(0, W-(1-R)), W)]   + [(0,W)]*2 + [(0,1)]*8 + [(0,10)]*7
    start = [0]*len(bounds)
    for l in range(len(bounds)):
        start[l] = np.random.uniform(bounds[l][0],bounds[l][1])

    eps = 1#e-3
    for _ in range(max_trials):
        result = opt.minimize(time, start, 
                bounds= bounds, tol=1e-6, 
                constraints=mycons, options={'maxiter':5000})
        
        adic= set_BJMM2even(*result.x)._asdict()
        if result.success and  abs( 2*adic['l2'] + adic['r1'] - adic['r0'] - adic['l0']) < eps: #2*x.l2 + x.r1 -x.r0, x.l0
            break
        else:
            for l in range(len(bounds)):
                start[l] = np.random.uniform(bounds[l][0],bounds[l][1])
    astuple = set_BJMM2even(*result.x)

    if verb:
        print("Validity: ", result.success,"|l0-l1|=",abs( 2*adic['l2'] + adic['r1'] - adic['r0'] - adic['l0']))
        print("q-Time: ", round_upwards_to_str(time(astuple)/scale/Q))
        print("2-Time: ", round_upwards_to_str(time(astuple)/scale))
        for t in adic:
            print(t, round_to_str(adic[t]) )
        print("Checking that the constraints are satisfied:")
        print(check_constraints(mycons, result.x))

    return time(astuple)/scale, adic, result.success



def make_BJMM2even_T_curve(q, z, R, prec, shift, onlyMax, avg, max_trials):
    Q = log(q, 2)
    Z = log(z, 2)
    
    E = get_E(q,z)

    if shift:
        shift = E[0]
        E = shift_E(E, shift)
        print("% shifted")

    if avg: print("%avg repr")
    else:   print("%max repr")

    if onlyMax: print("%only Max")
    else:   print("%all Diff")

    D = get_Diff(E, onlyMax)
    zd = len(D)

    nd = get_nd(E, True, avg)
    nb = get_nb(E, D, True, avg)
    print("% zd=",zd,"nd=",nd,"nb=",nb)
    
    print("%q=",q,"z=",z, "R=",round_to_str(R))
    for w in range(1,prec):
        W = w/prec
        if hbin(W) + W*Z < (1-R)*Q:
            if not shift or W >= 0.9:
                Cost2, params, success = optimize_BJMM2even(q, z, shift, onlyMax, avg, R, W, max_trials, verb=False)
                if success:
                    print("(",round_to_str(W),",",round_to_str(Cost2),") %", round_to_str(Cost2/Q))