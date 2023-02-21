"""
Simple script for calculating u_qm value of a fiber mode by solving
generalizezd characteristic equation in circular step index fiber
(infinite cladding).

author: jkkow
created: 2018-01-31T06:40:57
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, kv, jvp, kvp
from scipy.optimize import fsolve

V = 10
u = np.linspace(0, V, 10000)
w = np.sqrt(V*V-u*u)

q = 2
# n1/n2
n_ratio = 1.3
n_rs = n_ratio*n_ratio

############# Differential form of eigenvalue equation#
# TE = (u*w*w*jvp(q,u)/jv(q,u)+u*u*w*kvp(q,w)/kv(q,w))
# TM = (n_ratio*n_ratio*u*w*w*jvp(q,u)/jv(q,u)+u*u*w*kvp(q,w)/kv(q,w))
# LH = TE*TM
# RH = q*q*V*V*(n_ratio*n_ratio*w*w+u*u)


# for TE, TM, EH mode
TE_eh = u*w*w*jv(q+1, u)/jv(q, u) + u*u*w*kv(q+1, w)/kv(q, w) + q*V*V
TM_eh = n_rs*u*w*w*jv(q+1, u)/jv(q, u) + u*u*w*kv(q+1, w)/kv(q, w) + q*n_rs*u*u+w*w
LH_eh = np.sqrt(TE_eh)*np.sqrt(TM_eh)
RH_eh = q*V*np.sqrt(n_rs*u*u+w*w)

# for HE mode
TE_he = u*w*w*jv(q-1,u)/jv(q,u) - u*u*w*kv(q-1,w)/kv(q,w) - q*V*V
TM_he = n_rs*u*w*w*jv(q-1,u)/jv(q,u) - u*u*w*kv(q-1,w)/kv(q,w) - q*n_rs*u*u+w*w
LH_he = np.sqrt(TE_he)*np.sqrt(TM_he)
RH_he = -q*V*np.sqrt(n_rs*u*u+w*w)

LH_s = TE_eh*TM_eh
RH_s = q*q*V*V*(n_rs*u*u+w*w)


def characteristic(u):
    w = np.sqrt(V*V-u*u)
    TE_eh = u*w*w*jv(q+1,u)/jv(q,u) + u*u*w*kv(q+1,w)/kv(q,w) + q*V*V
    TM_eh = n_rs*u*w*w*jv(q+1,u)/jv(q,u) + u*u*w*kv(q+1,w)/kv(q,w) + q*n_rs*u*u+w*w
    LH_s = TE_eh*TM_eh
    RH_s = q*q*V*V*(n_rs*u*u+w*w)
    return LH_s - RH_s


def Brute_force_root_finder(f, start, end, steps):
    x = np.linspace(start, end, steps)
    y = f(x)
    roots = []
    for i in range(steps-1):
        if y[i]*y[i+1] < 0:
            root = x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])
            roots.append(root)
    return roots


if __name__ == "__main__":
    print("hello!")
    fig, axe = plt.subplots()
    # axe.plot(u, LH_weak, u, RH_weak)
    axe.plot(u,LH_eh , u, RH_eh)
    axe.set_ylim(-500,5000)
    fig2, axe2 = plt.subplots()
    axe2.plot(u, LH_s, u, RH_s)
    axe2.set_ylim(-2000, 1000000)

    roots = Brute_force_root_finder(characteristic, 1, V, 100000)
    print(roots)



    plt.show()
