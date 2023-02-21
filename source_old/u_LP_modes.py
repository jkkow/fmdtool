# -*- coding: utf-8 -*-
"""
Simple script for calculating u value of a fiber mode by solving
characteristic equation

author: jkkow
created: 2018-01-02T08:08:04
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jv, kv
from scipy.optimize import fsolve

V = 2.95262
u = np.linspace(0.001, V, 500)
w = np.sqrt(V * V - u * u)
l = 0

LH = u * jv(l - 1, u) / jv(l, u)
RH = -w * kv(l - 1, w) / kv(l, w)

fig, axe = plt.subplots()
axe.plot(u, LH, u, RH)
axe.plot(np.array(0))
axe.set_ylim(-15, 15)


def characteristic(u, V, l):
    w = np.sqrt(V * V - u * u)
    return u * jv(l - 1, u) / jv(l, u) + w * kv(l - 1, w) / kv(l, w)


def u_solve(u):
    return characteristic(u, V, l)


print("u%s1 = " % l, fsolve(u_solve, 1.5))

plt.show()
