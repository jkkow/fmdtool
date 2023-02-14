"""
This script returns a dictionary of u value in the given V-value.
For this, you should solve eigenvalue equation of LP mode first.
"""
from scipy.special import jv, kv, jn_zeros
import numpy as np


# The Eigenvalue equation to solve for the LP modes
def EEtoSolve(u, V, l):
    w = np.sqrt(V * V - u * u)
    return jv(l - 1, u) / jv(l, u) + (w / u) * kv(l - 1, w) / kv(l, w)


# Returns the list of diverging points
# (i.e. the u where J_l(u)=0) of the Eigenvalue equation
# in the given V and l value.
def DivPnts(V, l):
    div_pnts = []
    i = 1
    while jn_zeros(l, i)[-1] < V:
        div_pnts.append(jn_zeros(l, i)[-1])
        i += 1
    return div_pnts


# find maximum l value that has solution of the eigenvalue equation
def Maxl(V):
    l = 0  # start from l=0
    while DivPnts(V, l) != []:
        l += 1
    # if there's no diverging point in the eignevalue eq.,
    # check if the eigenvalue eq. at u=V larger than 0. If so, reduce 1 for l
    if jv(l - 1, V) / jv(l, V) > 0:
        l -= 1


def RootFind_startPnts(V):
    print()
    l_max = Maxl(V)
    rfsp = [[] for _ in range(l_max+1)]
    for l in range(l_max+1):
        print(l)
        rfsp[l] = [x-0.1 for x in DivPnts(V, l)]
        if jv(l - 1, V) / jv(l, V) < 0 and DivPnts(V, l) == []:
            rfsp[l].append(V/100)
        elif jv(l - 1, V) / jv(l, V) < 0:
            rfsp[l].append(DivPnts(V, l)[-1])

    return rfsp


def Sol_N(V):
    maxl = Maxl(V)
    sol_N = 0
