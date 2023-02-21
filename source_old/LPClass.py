"""
author: jkkow
created: 2018-05-06
updated: 2019-04-17
"""

import sys
import numpy as np
from scipy.special import jv, kv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation


class modeInspect(object):
    def __init__(self, V=10, a=10, u_dict=None):
        self.V = V
        self.a = a  # core radius
        if u_dict is not None:
            self.u_dict = u_dict
        else:
            self.u_dict = self.uCal(self.V, self.a)

        self.set_meshgrids(self.a)
        self.anim = []  # to handle multiple animation window

    def uCal(self, V, a):
        # Currently temporal fucntion. need to be finished afterward.
        udict = None
        assert (udict is not None), "Complete the fucntion 'uCal(V, a)'."
        return udict

    def set_meshgrids(self, a):
        # meshgrid for intensity plot
        i = np.linspace(-1.5 * a, 1.5 * a, 800)
        self.X, self.Y = np.meshgrid(i, i)
        self.r = np.sqrt(self.X**2 + self.Y**2)
        self.phi = np.arctan2(self.Y, self.X)
        # meshgrid for field vector plot
        i = np.linspace(-1.5 * a, 1.5 * a, 20)
        self.Xv, self.Yv = np.meshgrid(i, i)
        self.rv = np.sqrt(self.Xv**2 + self.Yv**2)
        self.phiv = np.arctan2(self.Yv, self.Xv)  # self.Yv must to come first

    # test functions for core region.
    # Rerturns 'True' for region, and  returns 'False' for non-interesting
    # region.
    def core(self, X, Y): return np.sqrt(X * X + Y * Y) <= self.a

    def clad(self, X, Y): return np.sqrt(X * X + Y * Y) > self.a

    def e_field_LP(self, l, m, r, phi, X, Y):
        """
        x-polarized field of LP mode.
        """
        u = self.u_dict[f'u{l}{m}']
        w = np.sqrt(self.V**2 - u**2)
        E_core = jv(l, u * r / self.a) * np.cos(l * phi)
        E_clad = (jv(l, u) / kv(l, w)) * \
            kv(l, w * r / self.a) * np.cos(l * phi)
        E_tot = self.core(X, Y) * E_core + self.clad(X, Y) * E_clad
        return E_tot

    def LP(self, l, m, ratio=1, pol=0, rot=0):
        if rot != 0:
            rad = rot * np.pi / 180
            phi = self.phi - rad
            phiv = self.phiv - rad
        else:
            phi = self.phi
            phiv = self.phiv

        E = self.e_field_LP(l, m, self.r, phi, self.X, self.Y)
        Ev = self.e_field_LP(l, m, self.rv, phiv, self.Xv, self.Yv)

        # normalization
        IntegralI = np.sum(np.abs(E * E))
        IntegralIv = np.sum(np.abs(Ev * Ev))
        # set total power to be 1
        Enorm = E * np.sqrt(ratio / IntegralI)
        Evnorm = Ev * np.sqrt(ratio / IntegralIv)

        # rotation of vector using vector rotation matrix. Here, Enorm and Evnorm
        # are actually x polarized field so we don't need to calculate full matrix
        # component.
        Ex = Enorm * np.cos(pol * np.pi / 180)
        Ey = Enorm * np.sin(pol * np.pi / 180)
        Exv = Evnorm * np.cos(pol * np.pi / 180)
        Eyv = Evnorm * np.sin(pol * np.pi / 180)

        return np.array([Ex, Ey, Exv, Eyv])

    def Gaussian_beam(self, ratio=1, pol=0, phase_diff=6.2832, beam_radi=None):
        """
        function Gaussian_beam() returns artificial Gaussian field array.
        Don't use as a reference for real Gaussian beam.

        <arguments>
        'phase_diff': sets phase shift differece between center and beam edge.
            the beam edge is determined by beam_radius.

        'beam_radi': sets beam radius. If it is not given, core radius is given.
        """

        if beam_radi is not None:
            waist = beam_radi
        else:
            waist = self.a

        r = np.sqrt(self.X**2 + self.Y**2)
        rv = np.sqrt(self.Xv**2 + self.Yv**2)
        C = phase_diff / (waist * waist)  # constant for radial phase shift

        # choose proper sign of radial phase term.
        Gaussian_E = np.exp(-(r / waist)**2) * np.exp(-1j * C * r * r)
        Gaussian_Ev = np.exp(-(rv / waist)**2) * np.exp(-1j * C * rv * rv)

        IntegralI = np.sum(np.abs(Gaussian_E**2))
        IntegralIv = np.sum(np.abs(Gaussian_Ev**2))

        Enorm = Gaussian_E * np.sqrt(ratio / IntegralI)
        Evnorm = Gaussian_Ev * np.sqrt(ratio / IntegralIv)

        Ex = Enorm * np.cos(pol * np.pi / 180)
        Ey = Enorm * np.sin(pol * np.pi / 180)
        Exv = Evnorm * np.cos(pol * np.pi / 180)
        Eyv = Evnorm * np.sin(pol * np.pi / 180)

        return np.array([Ex, Ey, Exv, Eyv])

    def plane_wave(self, ratio=1, pol=0):
        tEx = np.full_like(self.X, ratio / np.size(self.X))
        tExv = np.full_like(self.Xv, ratio / np.size(self.Xv))
        tEy = np.zeros_like(self.Y)
        tEyv = np.zeros_like(self.Yv)

        Ex = tEx * np.cos(pol * np.pi / 180)
        Ey = tEx * np.sin(pol * np.pi / 180)
        Exv = tExv * np.cos(pol * np.pi / 180)
        Eyv = tExv * np.sin(pol * np.pi / 180)

        return np.array([Ex, Ey, Exv, Eyv])

    def power(self, modeArray):
        Ex, Ey = modeArray[0], modeArray[1]
        P = np.sum(abs(Ex * Ex) + abs(Ey * Ey))
        return P

    def linear_polarizer(self, modeArray, deg):
        # upack mode information from mode array
        tEx, tEy = modeArray[0], modeArray[1]
        tExv, tEyv = modeArray[2], modeArray[3]

        Ex = (tEx * np.cos(deg * np.pi / 180.0) + tEy *
              np.sin(deg * np.pi / 180.0)) * np.cos(deg * np.pi / 180.0)
        Ey = (tEx * np.cos(deg * np.pi / 180.0) + tEy *
              np.sin(deg * np.pi / 180.0)) * np.sin(deg * np.pi / 180.0)
        Exv = (tExv * np.cos(deg * np.pi / 180.0) + tEyv *
               np.sin(deg * np.pi / 180.0)) * np.cos(deg * np.pi / 180.0)
        Eyv = (tExv * np.cos(deg * np.pi / 180.0) + tEyv *
               np.sin(deg * np.pi / 180.0)) * np.sin(deg * np.pi / 180.0)

        return np.array([Ex, Ey, Exv, Eyv])

    def birefrigence(self, modeArray, delphi, slow_axis='x'):
        # upack mode information from mode array
        tEx, tEy = modeArray[0], modeArray[1]
        tExv, tEyv = modeArray[2], modeArray[3]

        if slow_axis == 'x':
            Ex = tEx * np.exp(1j * delphi)
            Ey = tEy
            Exv = tExv * np.exp(1j * delphi)
            Eyv = tEyv
        elif slow_axis == 'y':
            Ex = tEx
            Ey = tEy * np.exp(1j * delphi)
            Exv = tExv
            Eyv = tEyv * np.exp(1j * delphi)
        else:
            raise Exception(
                "wrong 'slow_axis' argumet in function 'birefringence()'.")
        return np.array([Ex, Ey, Exv, Eyv])

    # gives phase shift in radian to mode array
    def delay(self, modeArray, shift):
        return modeArray * np.exp(1j * shift)

    def lobe_rot90(self, modeArray):
        # right agle rotation of mode.
        # upack mode information from mode array and change x, y values
        Ex = np.rot90(modeArray[0], 1)
        Ey = np.rot90(modeArray[1], 1)
        Exv = np.rot90(modeArray[2], 1)
        Eyv = np.rot90(modeArray[3], 1)

        return np.array([Ex, Ey, Exv, Eyv])

    def add_radial_phase(self, modeArray, phase_diff):
        """pahse_diff: radial phase difference between core center and boundary."""
        k = phase_diff / self.a
        r = np.sqrt(self.X**2 + self.Y**2)
        rv = np.sqrt(self.Xv**2 + self.Yv**2)

        Ex = np.exp(1j * k * r) * modeArray[0]
        Ey = np.exp(1j * k * r) * modeArray[1]
        Exv = np.exp(1j * k * rv) * modeArray[2]
        Eyv = np.exp(1j * k * rv) * modeArray[3]

        return np.array([Ex, Ey, Exv, Eyv])

    def mshow(self, modeArray, *args):
        # upack mode information from mode array
        Ex, Ey = modeArray[0], modeArray[1]
        Exv, Eyv = modeArray[2], modeArray[3]

        I = abs(Ex * Ex) + abs(Ey * Ey)  # intensity

        # set figure
        fig, ax = plt.subplots()
        ax.set_aspect("equal")

        plot_dict = {
            # without 'lambda', every function will be called when dictionary
            # called.
            'I': lambda: ax.pcolormesh(self.X, self.Y, I, cmap=cm.jet),
            'v': lambda: ax.quiver(self.Xv, self.Yv, np.real(Exv), np.real(Eyv), angles='xy', scale_units='xy', scale=0.3, pivot="middle", headwidth=5, headlength=3, headaxislength=2),
            'c': lambda: ax.contour(self.X, self.Y, I, levels=15, cmap=cm.afmhot),
            'A': lambda: ax.pcolormesh(self.X, self.Y, (np.real(Ex) + np.real(Ey)) / 2),
            'phase_x': lambda: ax.pcolormesh(self.X, self.Y, np.angle(Ex), cmap=cm.twilight_shifted),
            'phase_y': lambda: ax.pcolormesh(self.X, self.Y, np.angle(Ey), cmap=cm.twilight_shifted)
        }

        im = None  # to make it possible to plot only vector field
        for k in args:  # handling 'args' with 'plot_dict'
            if k == 'v':
                plot_dict.get(k, None)()
            else:
                im = plot_dict.get(k, None)()

        if im is not None:
            cb = plt.colorbar(im, orientation='vertical', spacing='uniform')
            # cb.set_clim(np.min(I), np.max(I))

        fig.tight_layout()
        plt.draw()

    def pmovie(self, modeArray):

        # upack mode information from mode array
        Ex, Ey = modeArray[0], modeArray[1]
        Exv, Eyv = modeArray[2], modeArray[3]

        I = abs(Ex * Ex) + abs(Ey * Ey)  # Intensity

        fig, ax = plt.subplots(figsize=(6, 6), facecolor='black')
        ax.set_aspect("equal")
        im = ax.pcolormesh(self.X, self.Y, I, cmap=cm.afmhot)
        ax.set_xticks([])
        ax.set_yticks([])

        # axins = inset_axes(ax,
        #                    width="5%",  # width = 10% of parent_bbox width
        #                    height="70%",  # height : 50%
        #                    loc=3,
        #                    bbox_to_anchor=(1.05, 0., 1, 1),
        #                    bbox_transform=ax.transAxes,
        #                    borderpad=0,
        #                    )
        #
        # cbar = fig.colorbar(im, cax=axins)
        # cbar.set_clim(0, 0.000015)
        # cbar.ax.tick_params(colors='white')
        # fig.tight_layout()

        qax = ax.quiver(self.Xv, self.Yv, np.real(Exv), np.real(Eyv), angles='xy', scale_units='xy', scale=0.3, pivot='tail', width=0.003,
                        headwidth=4, headlength=3, headaxislength=2)

        frame_num = 500
        delay = 10  # delay time btw frames in ms unit.

        def update_quiver(i, qax, Exv, Eyv):
            # time evolution of the x components of E field
            Exv = np.real(Exv * np.exp(-1j * i * 6 * np.pi / frame_num))
            # time evolution of the y components of E field
            Eyv = np.real(Eyv * np.exp(-1j * i * 6 * np.pi / frame_num))
            qax.set_UVC(Exv, Eyv)  # update U, V set
            return qax,

        # You have to put 'self' to the object variable of 'FuncAnimation()'
        # In case of multieple movies, list of 'self.anim' is configured.
        self.anim.append(
            animation.FuncAnimation(
                fig,
                update_quiver,
                frames=frame_num,
                fargs=(qax, Exv, Eyv),
                interval=delay,
                repeat_delay=delay,
                repeat=True,
                blit=True))

        plt.draw()


if __name__ == "__main__":
    # example: u values of V=10
    udict = {'u01': 2.18452, 'u02': 4.99665, 'u03': 7.76420, 'u11': 3.47699,
             'u12': 6.33103, 'u13': 9.04628, 'u21': 4.65441, 'u22': 7.56661, 'u31': 5.77402, 'u32': 8.72903, 'u41': 6.85500, 'u42': 9.81532}

    mi = modeInspect(u_dict=udict)

    lp11x = mi.LP(1, 1)
    mi.mshow(lp11x, 'phase_x')
    print(np.min(np.angle(lp11x[0])))
    print(np.angle(lp11x[2]))
    print(np.abs(lp11x[1]))
    plt.show()
