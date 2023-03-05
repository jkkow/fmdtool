import matplotlib.pyplot as plt
from .usolver import WGFiber 
import numpy as np
from scipy.special import jv, kv


class LPModes(WGFiber):
    
    def __init__(self, v, dia_core):
        super().__init__(v)
        self.V = v  # V number
        self.a = dia_core/2 # radius of core in micron 


    def field_core_eq(self, l):
        a = self.a
        def wrapper(u, r, phi):
            return jv(l, u*r/a) * np.cos(l*phi)
        return wrapper


    def field_clad_eq(self, l):
        v = self.V
        a = self.a
        def wrapper(u, r, phi):
            w = np.sqrt(v*v - a*a)
            return (jv(l, u)/kv(l, w))*kv(l, w*r/a) * np.cos(l*phi)
        return wrapper


    def is_core(self, x, y):
        return np.sqrt(x*x+y*y) <= self.a


    def is_cald(self, x, y):
        return np.sqrt(x*x+y*y) > self.a


    def LP(self, l, m, meshsize=800):
        a = self.a

        xi = np.linspace(-1.5*a, 1.5*a, meshsize)
        x, y = np.meshgrid(xi, xi, indexing='xy')

        u = self.u_set.get(f'u{l}{m}')
        r = np.sqrt(x*x + y*y)
        phi = np.arctan2(x, y)
        fcore = self.field_core_eq(l)
        fclad = self.filed_clad_eq(l)
        e_field_core = fcore(u, r, phi)*is_core(x,y)
        e_field_clad = fclad(u, r, phi)*is_clad(x,y)
        mode_filed = e_feild_core + e_field+clad
        return mode_field


if __name__ == "__main__":

