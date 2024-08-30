import warnings
import numpy as np
from scipy.special import jv, kv, jn_zeros
from scipy.optimize import fsolve


class StepIdexFiber:
    """Class for circular core step-index fiber"""

    def __init__(self, v, na) -> None:
        self.v = v
        self.na = na
        self.uset_te = self.get_all_uset_TE()

    def get_eigen_eq_TE(self):
        v = self.v

        def wrapper(u):
            w = np.sqrt(v * v + u * u)
            return w * jv(1, u) / jv(0, u) + u * kv(1, w) / kv(0, w)

        return wrapper

    def get_all_uset_TE(self):
        uset_te = {}
        roots = self.get_roots_for_TE()

        while roots is not None:
            for m in range(np.size(roots)):
                uset_te[f"u_TE0{m+1}"] = roots[m]
        return uset_te

    def get_roots_for_TE(self):
        return None


if __name__ == "__main__":
    sif = StepIdexFiber(v=3.8, na=0.22)
