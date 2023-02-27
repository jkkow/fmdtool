import warnings
import numpy as np
from scipy.special import jv, kv, jn_zeros
from scipy.optimize import fsolve


class WGFiber:

    def __init__(self, v=None):
        if v is not None:
            self.v = v


    def gen_eigen_eq(self, l):
        v = self.v
        def wrapper(u):
            w = np.sqrt(v*v - u*u)
            return u*jv(l-1, u)/jv(l, u) + w*kv(l-1, w)/kv(l, w)
        return wrapper


    @classmethod
    def left_side_eigen_eq(cls, l):
        def wrapper(u):
            return u*jv(l-1, u)/jv(l, u)
        return wrapper
        

    def find_max_jn_zeros(self, l):
        lo = 1
        while jn_zeros(l, lo)[-1] < self.v:
            lo += 1
            if lo > 100:
                raise ValueError("Too many points of Bessel zeros")
        maxn = lo
        return maxn-1


    def get_init_points_to_solve(self, l):
        v = self.v
        offset1 = 0.000001 # Don't change this value
        offset2 = 0.5
        lhs_eq_at_v = WGFiber.left_side_eigen_eq(l)(v)
        maxn = self.find_max_jn_zeros(l)
        if maxn == 0 :
            if lhs_eq_at_v > 0:
                return None
            else:
                init_points = np.array(v-offset1)
                return init_points
        else:
            init_points = jn_zeros(l, maxn) - offset2
            if lhs_eq_at_v > 0:
                return init_points
            else:
                return np.append(init_points, v-offset1)


    def get_roots_for_u(self, l):
        init_points = self.get_init_points_to_solve(l)
        if np.size(init_points) < 2 and init_points is None:
            return None
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            try:
                roots = fsolve(self.gen_eigen_eq(l), init_points)
            except RuntimeWarning as e:
                print(f"RuntimeWarning: {e}")
                print("The last solution is replaced by a value near V.")
                roots = fsolve(self.gen_eigen_eq(l), init_points[:-1])
                roots = np.append(roots, init_points[-1])
            finally:
                return roots


    def get_all_mode_set(self, l):
        pass


    @staticmethod
    def get_cuttoff_value(l, m):
        print(f"LP{l}{m}")
        cutoff_test_func = WGFiber.left_side_eigen_eq(l)(v)
    

        

if __name__ == "__main__":
    v = 7.02
    l = 0
    wgf = WGFiber(v)
    print(f"value at v of left-side = {wgf.left_side_eigen_eq(l)(wgf.v)} ")
    print(f"Init-points to solve LP{l}m = {wgf.get_init_points_to_solve(l)}")
    print(f"roots for LP{l}m mode = {wgf.get_roots_for_u(l)}")
    print(type(wgf.get_roots_for_u(l)))
    print(WGFiber.get_cuttoff_value(2,1))
   


