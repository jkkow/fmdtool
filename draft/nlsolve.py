from numpy import sqrt
from scipy.optimize import fsolve
from scipy.special import jv, kv

def eq_to_solve(V, u):
    w = sqrt(V*V-u*u)
    return 
