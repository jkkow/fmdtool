import warnings
import numpy as np
from scipy.special import jv, kv, jn_zeros
from scipy.optimize import fsolve


class WGFiber:
    """Class for Weakly Guiding Fiber"""

    def __init__(self, v):
        self.v = v
        self.uset = self.get_all_uset()

    def gen_eigen_eq(self, l):
        v = self.v

        def wrapper(u):
            w = np.sqrt(v * v - u * u)
            return u * jv(l - 1, u) / jv(l, u) + w * kv(l - 1, w) / kv(l, w)

        return wrapper

    @staticmethod
    def left_side_eigen_eq(l):
        """for some checks, plottings"""

        def wrapper(u):
            return u * jv(l - 1, u) / jv(l, u)

        return wrapper

    @staticmethod
    def right_side_eigen_eq(v, l):
        """for some checks, plottings"""

        def wrapper(u):
            w = np.sqrt(v * v - u * u)
            return -w * kv(l - 1, w) / kv(l, w)

        return wrapper

    def find_max_jn_zeros(self, l):
        """Find the number of points where a Bessel function
        jv(l, u) has the value zero in a range u is less than V"""
        nth = 1
        while jn_zeros(l, nth)[-1] < self.v:
            nth += 1
            if nth > 100:
                raise ValueError("Too many points of Bessel zeros")
        maxn = nth
        return maxn - 1

    def get_init_points_to_solve(self, l):
        v = self.v
        offset1 = 0.000001  # Don't change this value
        offset2 = 0.5
        lhs_eq_at_v = WGFiber.left_side_eigen_eq(l)(v)
        maxn = self.find_max_jn_zeros(l)
        if maxn == 0:
            if lhs_eq_at_v > 0:
                return None
            else:
                init_points = np.array(v - offset1)
                return init_points
        else:
            init_points = jn_zeros(l, maxn) - offset2
            if lhs_eq_at_v > 0:
                return init_points
            else:
                return np.append(init_points, v - offset1)

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

    def get_all_uset(self):
        uset = {}
        l = 0
        roots = self.get_roots_for_u(l)

        while roots is not None:
            for m in range(np.size(roots)):
                uset[f"u{l}{m+1}"] = roots[m]

            l += 1
            roots = self.get_roots_for_u(l)
        return uset

    def u_lm(self, l, m):
        if l < 0 or not isinstance(l, int):
            raise ValueError("'l' should be positive integers")
        if m < 1 or not isinstance(m, int):
            raise ValueError("'m' should be positive integers")

        try:
            return self.uset[f"u{l}{m}"]
        except KeyError as e:
            return None

    @staticmethod
    def get_cutoff_LP(l, m):
        if l < 0 or not isinstance(l, int):
            raise ValueError("'l' should be positive integers")

        if m < 1 or not isinstance(m, int):
            raise ValueError("'m' should be positive integers")

        if l == 0:
            try:
                return jn_zeros(l - 1, m)[-2]
            except IndexError as e:
                return 0.0
        else:
            return jn_zeros(l - 1, m)[-1]


if __name__ == "__main__":
    wgf = WGFiber(8.4174)
    print(wgf.get_all_uset())
