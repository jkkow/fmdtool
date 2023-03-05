import matplotlib.pyplot as plt
from usolver import WGFiber 
import numpy as np
from scipy.special import jv, kv


class LPModes(WGFiber):
    
    def __init__(self, v, dia_core):
        super().__init__(v)
        self.V = v  # V number
        self.a = dia_core/2 # radius of core in micron 
        self.mesh_size = 800


    def field_core_eq(self, l):
        a = self.a
        def wrapper(u, r, phi):
            return jv(l, u*r/a) * np.cos(l*phi)
        return wrapper


    def field_clad_eq(self, l):
        v = self.V
        a = self.a
        def wrapper(u, r, phi):
            w = np.sqrt(v*v - u*u)
            return (jv(l, u)/kv(l, w))*kv(l, w*r/a) * np.cos(l*phi)
        return wrapper


    def mask_core(self, x, y):
        return np.sqrt(x*x+y*y) <= self.a


    def mask_clad(self, x, y):
        return np.sqrt(x*x+y*y) > self.a


    def LP(self, l, m):
        a = self.a

        xi = np.linspace(-1.5*a, 1.5*a, self.mesh_size)
        x, y = np.meshgrid(xi, xi, indexing='xy')

        u = self.u_set.get(f'u{l}{m}')
        r = np.sqrt(x*x + y*y)
        phi = np.arctan2(x, y)

        fcore = self.field_core_eq(l)
        mcore = self.mask_core(x,y)
        e_field_core = fcore(u, r, phi)*mcore

        fclad = self.field_clad_eq(l)
        mclad = self.mask_clad(x,y)
        e_field_clad = fclad(u, r, phi)*mclad

        mode_field = e_field_core + e_field_clad
        return mode_field
