import numpy as np


class WGFiber:

    def __init__(self, V):
        self.V = V
        self.u = np.linspace(0.0000001, V-0.0001, int(V*1000))
    
    def gen_eigen_eq(V, l):
        pass

