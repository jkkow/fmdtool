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


    def left_side_eigen_eq(self, l):
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
        lhs_eq_at_v = self.left_side_eigen_eq(l)(v)
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
                roots = fsolve(self.gen_eigen_eq(l), init_points[:-1])
                roots = np.append(roots, init_points[-1])
                print("The last solution is replaced by a value near V.")
            finally:
                return roots


    def get_all_mode_set(self, l):
        pass


    @staticmethod
    def get_lp_cutoff(l, m):
        return jn_zeros(l-1, m)[-1]
        

if __name__ == "__main__":
    print(f"cutoff for LP01 = {WGFiber.get_cutoff_value('lp01')}")
    print(f"cutoff for LP11 = {WGFiber.get_cutoff_value('lp11')}")
    print(f"cutoff for LP21 = {WGFiber.get_cutoff_value('lp21')}")
    print(f"cutoff for LP02 = {WGFiber.get_cutoff_value('lp02')}")
