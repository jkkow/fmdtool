"""
This scripts defines class that can pursue the characteristics of LP modes
in step index circular core(infinite cladding) fiber.

author: jkkow
created: 2018-05-07T09:48:19
"""
import numpy as np
from scipy.special import jv, kv
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class ModePursuer(object):
    def __init__(self, V=8, a=17):
        self.V = V      # fiber V parameter
        self.a = a  # core radius
        self.u_dict = self.set_fiber_param()

        # set meshgrid for intensity plot
        i = np.linspace(-2 * a, 2 * a, 1000)
        self.X, self.Y = np.meshgrid(i, i)
        self.r = np.sqrt(self.X**2 + self.Y**2)
        self.phi = np.arctan2(self.Y, self.X)

        # set meshgrid for field vector plot
        i = np.linspace(-2 * a, 2 * a, 30)
        self.Xv, self.Yv = np.meshgrid(i, i)
        self.rv = np.sqrt(self.Xv**2 + self.Yv**2)
        self.phiv = np.arctan2(self.Yv, self.Xv)

    # Internal functions for the Class `PursueMode()`
    # internal class functions are named with lower case letters
    def set_fiber_param(self):
        u_dict = {'u01': 2.1246,
                  'u02': 4.8658,
                  'u03': 7.4529,
                  'u11': 3.943,
                  'u12': 6.1425,
                  'u21': 4.5383,
                  'u22': 7.2941,
                  'u31': 5.6215,
                  'u41': 6.6617}
        return u_dict

    # test function for core region.
    # This rerturns zero in the region outside core.
    def core(self, X, Y): return np.sqrt(X * X + Y * Y) <= self.a

    # test function for clad region. This returns zero in core region.
    def clad(self, X, Y): return np.sqrt(X * X + Y * Y) > self.a

    def e_field(self, l, m, u, w, r, phi, x, y):
        E_core = jv(l, u * r / self.a) * np.cos(l * phi) * self.core(x, y)
        E_clad = (jv(l, u) / kv(l, w)) * kv(l, w * r / self.a) * np.cos(l * phi) * self.clad(x, y)
        return E_core + E_clad

    def LP(self, l, m, power=1, pol=0, rot=0):
        """
        Function `LP` returns an array that consist of mode objects 'Ex', 'Ey', 'Exv', and 'Eyv'
        """
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

        # Getting intensity integral(power) for the normalizatoin
        IntegralE = sum(sum(abs(E * E)))  # for intensity plot
        IntegralEv = sum(sum(abs(Ev * Ev)))  # for field vector plot

        # Normalization of E filed for the intensity plot
        E_normal = np.sqrt(power) * E / np.sqrt(IntegralE)
        Ex = E_normal * np.cos(pol * np.pi / 180)
        Ey = E_normal * np.sin(pol * np.pi / 180)

        # Normalization of E filed for the filed vector plot
        Ev_normal = np.sqrt(power) * Ev / np.sqrt(IntegralEv)
        Exv = Ev_normal * np.cos(pol * np.pi / 180)
        Eyv = Ev_normal * np.sin(pol * np.pi / 180)

        description = "test temp"
        return np.array([Ex, Ey, Exv, Eyv, description])

    def ModePlot(self, modeObject, mode='I', field_vector='off'):
        Ex = modeObject[0]
        Ey = modeObject[1]
        Exv = modeObject[2]
        Eyv = modeObject[3]
        mode_descrpt = modeObject[4]

        Intensity = abs(Ex)**2 + abs(Ey)**2  # mesh for intensity plot
        Amplitude = np.real(Ex) + np.real(Ey)  # mesh for amplitude plot
        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(111, aspect="equal")
        # ax.set_title(mode_descrpt, size=25)

        if field_vector == "off" and mode == "I":
            im = ax.pcolormesh(self.X, self.Y, Intensity, cmap=cm.afmhot)
        elif field_vector == "off" and mode == "A":
            im = ax.pcolormesh(self.X, self.Y, Amplitude, cmap=cm.jet)
        elif field_vector == "on" and mode == "I":
            # caution! 'quiver' function should be located after 'pcolormesh'.
            # Or it will be covered by result of 'pcolormesh' and no show.
            im = ax.pcolormesh(self.X, self.Y, Intensity, cmap=cm.afmhot)
            ax.quiver(self.Xv, self.Yv, Exv, Eyv, scale=4, pivot="middle")
        elif field_vector == "on" and mode == "A":
            im = ax.pcolormesh(self.X, self.Y, Amplitude, cmap=cm.jet)
            ax.quiver(self.Xv, self.Yv, Exv, Eyv, scale=4, pivot="middle")

        plt.colorbar(im)
        plt.draw()

    def PhaseDelay(self, modeObject, shift):
        pass

    def ModeMovie(self, modeObject):
        pass

    def ModeSum(self, modeObjects):
        pass

    def Polarizer(self, modeObject, pol_angle='x'):
        pass

    def ModeRotate(self, modeObject, rot_angle=0):
        pass

    def PossibleModes(self, printout='off'):
        pass


if __name__ == "__main__":
    mp = ModePursuer(V=8, a=17)

    LP01 = mp.LP(0, 1, power=0.5, pol=0, rot=0)
    # LP11 = mp.LP(1, 1)
    # LP21 = mp.LP(2, 1)
    # LP02 = mp.LP(0, 2)
    # LP21 = mp.LP(2, 1, power=1, pol=0, rot=0)
    LP02 = mp.LP(0, 2, power=1, pol=0, rot=0)
    LP11 = mp.LP(1,1)

    # LP31 = mp.LP(3, 1, power=0.1, pol=0, rot=0)
    # LP12 = mp.LP(1, 2, power=1, pol=0, rot=0)
    # mp.ModePlot(LP21, mode='A', field_vector='on')
    mp.ModePlot(LP11, mode='A')
    mp.ModePlot(LP01 + LP02, mode='I')
    mp.ModePlot(LP02, mode='I')

    # ms = LP21 + LP02
    # mp.ModePlot(ms, mode='I', field_vector='off')

    plt.show()

    # PossModeNum, ModesSet = mp.PossibleModes()
    # print(PossModeNum)
    # print(ModeSet)
    # mp.PssibleModes(printout='on')
    #
    # LP11 = mp.LP(1, 1, power=1, pol=0, rot=0)
    # LP21 = mp.LP(2, 1, power=1, pol=0, rot=0)
    # LP02 = mp.LP(0, 2, power=0.1, pol=0, rot=0)
    #
    # modeTotal = mp.ModeSum(LP11 + LP21)
    # mp.PhaseDelay(LP02, np.pi / 4)
    # modeTotal = mp.ModeSum(modeTotal + LP02)
    # mp.History(modeTotal)
    #
    # mp.ModePlot(LP11, mode='A', field_vector='on')
    # mp.ModePlot(LP11 + LP21)
    # shift = np.pi / 4
    # mp.ModePlot(LP11 + 1j * shift * LP21)
    #
    # mp.ModeMovie(LP21)
    # mp.ModeMovie(LP21 + 1j * shift * LP11)
    # mp.ModeMovie(modeTotal)
