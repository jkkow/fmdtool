from .usolver import WGFiber
import numpy as np
from scipy.special import jv, kv


class LPModes(WGFiber):
    def __init__(self, v, d):
        super().__init__(v)
        self.v = v  # V number
        self.a = d/ 2  # radius of core in micron
        self.mesh_size = 300
        self.xi = np.linspace(-2 * self.a, 2 * self.a, self.mesh_size)

    def field_core_eq(self, l, u):
        def wrapper(r, phi):
            return jv(l, u * r / self.a) * np.cos(l * phi)
        return wrapper

    def field_clad_eq(self, l, u):
        def wrapper(r, phi):
            w = np.sqrt(self.v * self.v - u * u)
            return (jv(l, u) / kv(l, w)) * kv(l, w * r / self.a) * np.cos(l * phi)

        return wrapper

    def mask_core(self, x, y):
        return np.sqrt(x * x + y * y) <= self.a

    def mask_clad(self, x, y):
        return np.sqrt(x * x + y * y) > self.a

    def gen_mode_LP(self, l, m, rot=0, jones=(1, 0)):
        # Returns complex wave function
        u = self.u_lm(l, m)  # Method of 'WGFiber'

        jones_vector = np.array(jones)
        norm_jones = np.sqrt(np.sum(jones_vector.real**2 + jones_vector.imag**2))
        if not np.isclose(norm_jones, 1.0):
            raise ValueError("Jones Vector should be normalized.")

        if u is not None:
            x, y = np.meshgrid(self.xi, self.xi, indexing="xy")
            r = np.sqrt(x * x + y * y)
            phi = np.arctan2(x, y) - rot*np.pi/180 # unit of rot is degree

            fcore = self.field_core_eq(l, u)
            mcore = self.mask_core(x, y)  # returns 1 for core region, 0 for clad region.
            e_field_core = fcore(r, phi) * mcore

            fclad = self.field_clad_eq(l, u)
            mclad = self.mask_clad(x, y)  # returns 0 for core region, 1 for clad region.
            e_field_clad = fclad(r, phi) * mclad

            field = e_field_core + e_field_clad
            norm_field = field / np.sqrt(np.sum(field.real**2 + field.imag**2))
            Ex = jones_vector[0]*norm_field
            Ey = jones_vector[1]*norm_field
            return np.array([Ex, Ey])
        else:
            raise ValueError(f"LP{l}{m} mode doesn't exist at V={self.v}.")

    def get_power(self, mobj):
        return np.sum(mobj.real**2 + mobj.imag**2)
