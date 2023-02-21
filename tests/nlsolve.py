from numpy import sqrt
from scipy.optimize import fsolve
from scipy.special import jv, kv
import matplotlib.pyplot as plt

def eq_to_solve(V, l):
    def func(u):
        w = sqrt(V*V - u*u)
        return jv(l - 1, u) / jv(l, u) + (w / u) * kv(l - 1, w) / kv(l, w)
    return func


if __name__ == "__main__":
    root = fsolve(eq_to_solve(V=8, l=0), [2.1, 5.2, 7.0])
    print(root)
