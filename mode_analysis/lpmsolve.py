import numpy as np
from scipy.optimize import fsolve
from scipy.special import jv, kv


class LPModes:

    def __init__(self, V):
        self.V = V


    def get_field_ex_core(V, u):
        w = np.sqrt(V*V - u*u)


class UlmSolver:

    def __init__(self, V):
        self.V = V
        self.u = np.linspace(0, V, 

    def get_eigen_eq(self, V, u, l):
        w = np.sqrt(V * V - u * u)
        return jv(l - 1, u) / jv(l, u) + (w / u) * kv(l - 1, w) / kv(l, w)

