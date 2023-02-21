"""
This LP mode class automatically generates new LP mode instance when pre-made mode instances are added so that pre-made modes can remain unaffected.
This approace is very intuitive to add an effect to eigen mode instance, but it takes long time to operate and is difficult to manupulate the filed of the mode instance

author: jkkow
created: 2018-05-06
"""


import numpy as np
from scipy.special import jv, kv
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class LP_mode_inspector(object):
    def __init__(self, V=8, a=17):
        self.V = V
        self.a = a  # core radius
        # dictionary of u value for each LP mode. Need to be calculated prior with corresponding V.
        self.u_dict = {'u01': 2.1246, 'u02': 4.8658, 'u03': 7.4529, 'u11': 3.943,
                       'u12': 6.1425, 'u21': 4.5383, 'u22': 7.2941, 'u31': 5.6215, 'u41': 6.6617}

        # meshgrid for intensity plot
        i = np.linspace(-2 * a, 2 * a, 500)
        self.X, self.Y = np.meshgrid(i, i)
        self.r = np.sqrt(self.X**2 + self.Y**2)
        self.phi = np.arctan2(self.Y, self.X)
        # meshgrid for field vector plot
        i = np.linspace(-2 * a, 2 * a, 30)
        self.Xv, self.Yv = np.meshgrid(i, i)
        self.rv = np.sqrt(self.Xv**2 + self.Yv**2)
        self.phiv = np.arctan2(self.Yv, self.Xv)

        # set 'Ex' has same dimnesion with 'self.r' and value of zero.
        self.Ex = 0 * self.r
        self.Ey = 0 * self.r
        self.Exv = 0 * self.rv
        self.Eyv = 0 * self.rv
        self.I = 0 * self.r
        self.A = 0 * self.rv
        self.description = ""

    # test function for core region.
    # This rerturns 'True' for inside core, 'False' for outsiede.
    def core(self, X, Y): return np.sqrt(X * X + Y * Y) <= self.a

    # test function for clad region.
    # This rerturns 'True' for ourside core, 'False' for inside.
    def clad(self, X, Y): return np.sqrt(X * X + Y * Y) > self.a

    # filed equation for LP modes in circular core step index fiber
    def e_field(self, l, m, u, w, r, phi, x, y):
        E_core = jv(l, u * r / self.a) * np.cos(l * phi) * self.core(x, y)
        E_clad = (jv(l, u) / kv(l, w)) * kv(l, w * r / self.a) * np.cos(l * phi) * self.clad(x, y)
        return E_core + E_clad

    def LP(self, l, m, ratio=1, pol=0, rot=0):
        u = self.u_dict["u{}{}".format(l, m)]
        w = np.sqrt(self.V**2 - u**2)

        if rot != 0:
            rad = rot * np.pi / 180
            phi = self.phi - rad
            phiv = self.phiv - rad
            E = self.e_field(l, m, u, w, self.r, phi, self.X, self.Y)
            Ev = self.e_field(l, m, u, w, self.rv, phiv, self.Xv, self.Yv)
        else:
            E = self.e_field(l, m, u, w, self.r, self.phi, self.X, self.Y)
            Ev = self.e_field(l, m, u, w, self.rv, self.phiv, self.Xv, self.Yv)

        IntegralE = sum(sum(abs(E * E)))
        # normailzed E field to genrate intensity of 'ratio'
        E_normal = np.sqrt(ratio) * E / np.sqrt(IntegralE)
        Ex = E_normal * np.cos(pol * np.pi / 180)
        Ey = E_normal * np.sin(pol * np.pi / 180)

        # E filed for the meshgrid of polaization plot
        IntegralEv = sum(sum(abs(Ev * Ev)))
        Ev_normal = np.sqrt(ratio) * Ev / np.sqrt(IntegralEv)
        Exv = Ev_normal * np.cos(pol * np.pi / 180)
        Eyv = Ev_normal * np.sin(pol * np.pi / 180)

        LPinstance = LP_mode_inspector(V=self.V)
        LPinstance.Ex = Ex
        LPinstance.Ey = Ey
        LPinstance.Exv = Exv
        LPinstance.Eyv = Eyv
        LPinstance.description = "LP$_{}$$_{}$".format(l, m)
        return LPinstance

    def Delay(self, shift):
        self.Ex = self.Ex * np.exp(1j * shift)
        self.Ey = self.Ey * np.exp(1j * shift)
        self.Exv = self.Exv * np.exp(1j * shift)
        self.Eyv = self.Eyv * np.exp(1j * shift)
        self.description += r"$^\mathtt{%.2f\pi}$" % float(shift / np.pi)

    def __add__(self, other):
        Ex = self.Ex + other.Ex
        Ey = self.Ey + other.Ey
        Exv = self.Exv + other.Exv
        Eyv = self.Eyv + other.Eyv
        newInstance = LP_mode_inspector(V=self.V)
        newInstance.Ex = Ex
        newInstance.Ey = Ey
        newInstance.Exv = Exv
        newInstance.Eyv = Eyv
        newInstance.description = self.description + "+ " + other.description
        return newInstance

    def ModePlot(self, mode="I", vector="off"):
        self.I = abs(self.Ex)**2 + abs(self.Ey)**2  # mesh for intensity plot
        self.A = np.real(self.Ex) + np.real(self.Ey)  # mesh for amplitude plot
        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(111, aspect="equal")
        ax.set_title(self.description, size=25)

        if vector == "off" and mode == "I":
            im = ax.pcolormesh(self.X, self.Y, self.I, cmap=cm.afmhot)
        elif vector == "off" and mode == "A":
            im = ax.pcolormesh(self.X, self.Y, self.A, cmap=cm.jet)
        elif vector == "on" and mode == "I":
            # caution! 'quiver' function should be located after 'pcolormesh'.
            # Or it will be covered by result of 'pcolormesh' and no show.
            im = ax.pcolormesh(self.X, self.Y, self.I, cmap=cm.afmhot)
            ax.quiver(self.Xv, self.Yv, self.Exv, self.Eyv, scale=4, pivot="middle")
        elif vector == "on" and mode == "A":
            im = ax.pcolormesh(self.X, self.Y, self.A, cmap=cm.jet)
            ax.quiver(self.Xv, self.Yv, self.Exv, self.Eyv, scale=4, pivot="middle")

        plt.colorbar(im)
        plt.draw()


if __name__ == '__main__':
    pm = LP_mode_inspector()

    LP11 = pm.LP(1, 1, rot=45)

    LP11.ModePlot(mode='A', vector='on')
    plt.show()
