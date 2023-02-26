import numpy as np
from scipy.special import jv, kv, jn_zeros
from scipy.optimize import fsolve


class WGFiber:


    def __init__(self, v):
        self.v = v

    def gen_eigen_eq(self, l):
        v = self.v
        def wrapper(u):
            w = np.sqrt(v*v - u*u)
            return u*jv(l-1, u)/jv(l, u) + w*kv(l-1, w)/kv(l, w)
        return wrapper


    def left_side_eigen_eq(self, l):
        def wrapper(u):
            return u*jv(l-1, u)/jv(l, u)
        return wrapper
        

    def find_max_jn_zeros(self, l):
        start = 1
        while jn_zeros(l, start)[-1] < self.v:
            start += 1
        max_num = start
        return max_num-1


    def get_diverging_points(self, l):
        maxnum = self.find_max_jn_zeros(l)
        return jn_zeros(l, maxnum)
        

    def usolve(self, l):
        eigen_eq = self.gen_eigen_eq(l)
        initp = self.get_init_points_to_solve(l)
        roots = fsolve(eigen_eq, initp)
        dict_u = {}
        for m in len(roots):
            dict_u[f"ul{m}"] = roots[m]
        return dict_u

            
if __name__ == "__main__":
    wgf = WGFiber(5)
