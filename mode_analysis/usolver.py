import numpy as np
from scipy.special import jv, kv


class WGFiber:

    def __init__(self, V=None):
        if V is not None:
            self.V = V

    def gen_eigen_eq(self, V, l):
        def wrapper(u):
            w = np.sqrt(V*V - u*u)
            return jv(l - 1, u) / jv(l, u) + (w / u) * kv(l - 1, w) / kv(l, w)
        return wrapper
