"""
functions used in our optimizations.
"""

from math import*
import scipy.optimize as opt


def check_constraints(constraints, solution) : 
    return [ (constraint['type'], constraint['fun'](solution)) for constraint in constraints ]


def wrap(f,g) :
    def inner(x):
        return f(g(*x))
    return inner



#=================================
# general definitions
#=================================

def xlx(x):
    if x<=0: return - 100*x
    return x*log(x, 2)


def hbin(x):
    """
    binary entropy function of x
    """
    if 0<x<1: return -x*log(x, 2) - (1-x)*log(1-x, 2)
    return 0

def trinom(a,b):
    """
    exponential coeff of trinomial
    """
    return - xlx(a) - xlx(b) - xlx(1-a-b)

def quad(a, b, c):
    return -xlx(a) - xlx(b) - xlx(c) - xlx(1-a-b-c)
#=================================
# macros used by stern
#=================================

def Niter(R,W,L,P):
    """
    determine the number of iterations as P_succ^-1
    R: rate
    W: overall weight
    L: red small instance
    P: small instance
    """
    return -xlx(W)-xlx(1-W) -xlx(1-R-L) + xlx(W-P) + xlx(1-R-L-W+P) - xlx(R+L) + xlx(P) + xlx(R+L-P)


def Lbase(R,L,P,Z):
    """
    base list for stern with wt p in small instance:
    b_i = ((k+l)/2, p/2) * z^(p/2)
    R: rate
    W: overall weight
    L: red small instance
    P: small instance
    """
    return (R+L)/2 * hbin(P/(R+L)) + Z*P/2



def round_to_str(t):
    """
    Rounds the value 't' to a string with 4 digit precision (adding trailing zeroes
    to emphasize precision).
    """
    s = str(round(t,4))
    # must be 6 digits
    return (s + "0" * (5 + s.find(".") -len(s)))
    #return t   

