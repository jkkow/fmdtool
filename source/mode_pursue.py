"""
Created on Sat Aug 15 21:43:43 2015

modules that pursues effects of mode manupulations.
Conducting various calculations with modes in weakly guiding step index fiber, and visualize the effect.
@author: JKwonKo

need to be updated for amplitude normalization
"""
from pylab import *
from scipy.special import jv, kv
from matplotlib import animation

a = 10   # micrometer
V = 5
u01 = 1.9737
u02 = 4.3490
u11 = 3.1163
u21 = 4.1204

x1 = linspace(-2 * a, 2 * a, 30)
y1 = linspace(-2 * a, 2 * a, 30)
Xf, Yf = meshgrid(x1, y1)
x2 = linspace(-2 * a, 2 * a, 600)
y2 = linspace(-2 * a, 2 * a, 600)
X, Y = meshgrid(x2, y2)


def core(X, Y): return sqrt(X * X + Y * Y) <= a


def clad(X, Y): return sqrt(X * X + Y * Y) >= a


def phi(X, Y):
    def qdrt13(X, Y): return X * Y > 0

    def qdrt24(X, Y): return X * Y < 0

    def qdrt14(X, Y): return Y * Y * X > 0

    def qdrt23(X, Y): return Y * Y * X < 0
    qdrt14 = qdrt14(X, Y)
    qdrt2 = qdrt24(X, Y) * qdrt23(X, Y)
    qdrt3 = qdrt13(X, Y) * qdrt23(X, Y)
    phi14 = arctan(Y / X) * qdrt14
    phi2 = (arctan(Y / X) + pi) * qdrt2
    phi3 = (arctan(Y / X) - pi) * qdrt3
    phi = phi14 + phi2 + phi3
    return phi


core_f = core(Xf, Yf)
clad_f = clad(Xf, Yf)
phi_f = phi(Xf, Yf)
r_f = sqrt(Xf * Xf + Yf * Yf)

core = core(X, Y)
clad = clad(X, Y)
phi = phi(X, Y)
r = sqrt(X * X + Y * Y)


def u(l, m):
    if l == 0 and m == 1:
        return u01
    elif l == 0 and m == 2:
        return u02
    elif l == 1 and m == 1:
        return u11
    elif l == 1 and m == 2:
        return u12
    elif l == 1 and m == 3:
        return u13
    elif l == 2 and m == 1:
        return u21
    elif l == 3 and m == 1:
        return u31


def eq_he(q, u, w, r, phi, core, clad):
    Er_core = -jv(q - 1, u * r / a) * sin(q * phi) * core
    Er_clad = -(jv(q - 1, u) * kv(q - 1, w * r / a) * sin(q * phi) / kv(q - 1, w)) * clad
    Ephi_core = -jv(q - 1, u * r / a) * cos(q * phi) * core
    Ephi_clad = -(jv(q - 1, u) * kv(q - 1, w * r / a) * cos(q * phi) / kv(q - 1, w)) * clad

    Er = Er_core + Er_clad
    Ephi = Ephi_core + Ephi_clad
    Ex = Er * cos(phi) - Ephi * sin(phi)
    Ey = Er * sin(phi) + Ephi * cos(phi)
    return Er, Ephi, Ex, Ey


def eq_eh(q, u, w, r, phi, core, clad):
    Er_core = jv(q + 1, u * r / a) * sin(q * phi) * core
    Er_clad = (jv(q + 1, u) * kv(q + 1, w * r / a) * sin(q * phi) / kv(q + 1, w)) * clad
    Ephi_core = -jv(q + 1, u * r / a) * cos(q * phi) * core
    Ephi_clad = -(jv(q + 1, u) * kv(q + 1, w * r / a) * cos(q * phi) / kv(q + 1, w)) * clad

    Er = Er_core + Er_clad
    Ephi = Ephi_core + Ephi_clad
    Ex = Er * cos(phi) - Ephi * sin(phi)
    Ey = Er * sin(phi) + Ephi * cos(phi)
    return Er, Ephi, Ex, Ey


def eq_te(q, u, w, r, phi, core, clad):
    Er_core = 0
    Er_clad = 0
    Ephi_core = -jv(1, u * r / a) * core
    Ephi_clad = -(jv(1, u) * kv(1, w * r / a) / kv(1, w)) * clad

    Er = Er_core + Er_clad
    Ephi = Ephi_core + Ephi_clad
    Ex = Er * cos(phi) - Ephi * sin(phi)
    Ey = Er * sin(phi) + Ephi * cos(phi)
    return Er, Ephi, Ex, Ey


def eq_tm(q, u, w, r, phi, core, clad):
    Er_core = jv(1, u * r / a) * core
    Er_clad = (jv(1, u) * kv(1, w * r / a) / kv(1, w)) * clad
    Ephi_core = 0
    Ephi_clad = 0

    Er = Er_core + Er_clad
    Ephi = Ephi_core + Ephi_clad
    Ex = Er * cos(phi) - Ephi * sin(phi)
    Ey = Er * sin(phi) + Ephi * cos(phi)
    return Er, Ephi, Ex, Ey


class pursue:

    fEr = Xf * 0
    fEphi = Xf * 0
    fEx = Xf * 0
    fEy = Xf * 0

    Er = X * 0
    Ephi = X * 0
    Ex = X * 0
    Ey = X * 0

    history = ""

    def he(self, q, m):
        self.q = q
        self.m = m
        self.l = q - 1
        self.u = u(self.l, self.m)
        self.w = sqrt(V**2 - self.u**2)
        self.fEr, self.fEphi, self.fEx, self.fEy = eq_he(
            self.q, self.u, self.w, r_f, phi_f, core_f, clad_f)
        self.Er, self.Ephi, self.Ex, self.Ey = eq_he(self.q, self.u, self.w, r, phi, core, clad)
        self.history += " >> HE%d%d" % (self.q, self.m)

    def eh(self, q, m):
        self.q = q
        self.m = m
        self.l = q + 1
        self.u = u(self.l, self.m)
        self.w = sqrt(V**2 - self.u**2)
        self.fEr, self.fEphi, self.fEx, self.fEy = eq_eh(
            self.q, self.u, self.w, r_f, phi_f, core_f, clad_f)
        self.Er, self.Ephi, self.Ex, self.Ey = eq_eh(self.q, self.u, self.w, r, phi, core, clad)
        self.history += " >> EH%d%d" % (self.q, self.m)

    def te(self, q, m):
        if q != 0:
            self.q = 0
            print("q should be zero for TE mode. It'is automatically changed to 0.")
        self.q = q
        self.m = m
        self.l = 1
        self.u = u(self.l, self.m)
        self.w = sqrt(V**2 - self.u**2)
        self.fEr, self.fEphi, self.fEx, self.fEy = eq_te(
            self.q, self.u, self.w, r_f, phi_f, core_f, clad_f)
        self.Er, self.Ephi, self.Ex, self.Ey = eq_te(self.q, self.u, self.w, r, phi, core, clad)
        self.history += " >> TE%d%d" % (self.q, self.m)

    def tm(self, q, m):
        if q != 0:
            self.q = 0
            print("q should be zero for TM mode. It'is automatically changed to 0.")
        self.q = q
        self.m = m
        self.l = 1
        self.u = u(self.l, self.m)
        self.w = sqrt(V**2 - self.u**2)
        self.fEr, self.fEphi, self.fEx, self.fEy = eq_tm(
            self.q, self.u, self.w, r_f, phi_f, core_f, clad_f)
        self.Er, self.Ephi, self.Ex, self.Ey = eq_tm(self.q, self.u, self.w, r, phi, core, clad)
        self.history += " >> TM%d%d" % (self.q, self.m)

    def delay(self, shift):
        self.fEr = self.fEr * exp(1j * shift)
        self.fEphi = self.fEphi * exp(1j * shift)
        self.fEx = self.fEx * exp(1j * shift)
        self.fEy = self.fEy * exp(1j * shift)
        self.Er = self.Er * exp(1j * shift)
        self.Ephi = self.Ephi * exp(1j * shift)
        self.Ex = self.Ex * exp(1j * shift)
        self.Ey = self.Ey * exp(1j * shift)
        self.history += " >> delay(%.3f)" % shift

    def polarizer(self, degr, justcheck='on'):
        self.tem1f = (self.fEx * cos(degr * pi / 180.0) + self.fEy *
                      sin(degr * pi / 180.0)) * cos(degr * pi / 180.0)
        self.tem2f = (self.fEx * cos(degr * pi / 180.0) + self.fEy *
                      sin(degr * pi / 180.0)) * sin(degr * pi / 180.0)
        self.tem1 = (self.Ex * cos(degr * pi / 180.0) + self.Ey *
                     sin(degr * pi / 180.0)) * cos(degr * pi / 180.0)
        self.tem2 = (self.Ex * cos(degr * pi / 180.0) + self.Ey *
                     sin(degr * pi / 180.0)) * sin(degr * pi / 180.0)

        if justcheck == 'on':
            self.I = abs(self.tem1)**2 + abs(self.tem2)**2
            fig = figure()
            im = pcolormesh(X, Y, self.I, cmap=cm.afmhot)
            colorbar(im)
            quiver(Xf, Yf, real(self.tem1f), real(self.tem2f), pivot='middle', scale=20)
        else:
            self.fEx = self.tem1f
            self.fEy = self.tem2f
            self.fEr = self.fEx * cos(phi_f) + self.fEy * sin(phi_f)
            self.fEphi = -self.fEx * sin(phi_f) + self.fEy * cos(phi_f)
            self.Ex = self.tem1
            self.Ey = self.tem2
            self.Er = self.Ex * cos(phi) + self.Ey * sin(phi)
            self.Ephi = -self.Ex * sin(phi) + self.Ey * cos(phi)
            self.history += " >> polarizer(%.2f)" % degr

    def orthogonal(self):
        self.fEr, self.fEphi = -self.fEphi, self.fEr
        self.fEx = self.fEr * cos(phi_f) - self.fEphi * sin(phi_f)
        self.fEy = self.fEr * sin(phi_f) + self.fEphi * cos(phi_f)
        self.Er, self.Ephi = -self.Ephi, self.Er
        self.Ex = self.Er * cos(phi) - self.Ephi * sin(phi)
        self.Ey = self.Er * sin(phi) + self.Ephi * cos(phi)
        self.history += " >> orthogoal()"

    def __add__(self, other):
        self.fEr = self.fEr + other.fEr
        self.fEphi = self.fEphi + other.fEphi
        self.fEx = self.fEx + other.fEx
        self.fEy = self.fEy + other.fEy
        self.Er = self.Er + other.Er
        self.Ephi = self.Ephi + other.Ephi
        self.Ex = self.Ex + other.Ex
        self.Ey = self.Ey + other.Ey
        self.history += " >> add [%s]" % other.history

    def __sub__(self, other):
        self.fEr = self.fEr - other.fEr
        self.fEphi = self.fEphi - other.fEphi
        self.fEx = self.fEx - other.fEx
        self.fEy = self.fEy - other.fEy
        self.Er = self.Er - other.Er
        self.Ephi = self.Ephi - other.Ephi
        self.Ex = self.Ex - other.Ex
        self.Ey = self.Ey - other.Ey
        self.history += " >> substract [%s]" % other.history

    def birfg(self, delphi):
        self.fEx = self.fEx * exp(1j * delphi)
        self.fEr = self.fEx * cos(phi_f) + self.fEy * sin(phi_f)
        self.fEphi = -self.fEx * sin(phi_f) + self.fEy * cos(phi_f)
        self.Ex = self.Ex * exp(1j * delphi)
        self.Er = self.Ex * cos(phi) + self.Ey * sin(phi)
        self.Ephi = -self.Ex * sin(phi) + self.Ey * cos(phi)
        self.history += " >> birefringence(%.2f)" % delphi

    def draw(self, field='on', inten='on', amp='off'):

        self.I = abs(self.Ex)**2 + abs(self.Ey)**2
        self.A = real(self.Ex) + real(self.Ey)
        fig = figure()
        if field == 'on' and inten == 'on' and amp == 'off':
            im = pcolormesh(X, Y, self.I, cmap=cm.afmhot)
            colorbar(im)
            quiver(Xf, Yf, real(self.fEx), real(self.fEy), scale=20, pivot='middle')
            # fig.savefig('%s' % self.history, format='jpg', dpi=600)

        elif field == 'on' and inten == 'off' and amp == 'off':
            imshow(Xf * 0, cmap=cm.bone_r, aspect='equal', extent=[-2 * a, 2 * a, -2 * a, 2 * a])
            quiver(Xf, Yf, real(self.fEx), real(self.fEy), scale=20, pivot='middle')
        elif field == 'off' and inten == 'on' and amp == 'off':
            im = pcolormesh(X, Y, self.I, cmap=cm.afmhot)
            colorbar(im)
        elif field == 'on' and inten == 'off' and amp == 'on':
            im = pcolormesh(X, Y, self.A)
            colorbar(im)
            quiver(Xf, Yf, self.fEx, self.fEy, scale=20, pivot='middle')
        elif amp == 'on':
            im = pcolormesh(X, Y, self.A)
            colorbar(im)
        else:
            print("\nTry correct arguments.")

        draw()

    def movie(self):

        self.I = abs(self.Ex)**2 + abs(self.Ey)**2
        self.A = real(self.Ex) + real(self.Ey)

        fig = figure()
        im = pcolormesh(X, Y, self.I, cmap=cm.afmhot)
        colorbar(im)
        Q = quiver(Xf, Yf, Xf * 0, Yf * 0, scale=20)

        def update_quiver(num, Q, U, V):
            U = real(U * exp(-1j * num * 0.2))  # time evolution of the x components of vectors
            V = real(V * exp(-1j * num * 0.2))  # time evolution of the y components of vectors
            Q.set_UVC(U, V)  # update U,V set
            return Q,

        self.anima = animation.FuncAnimation(fig, update_quiver, fargs=(
            Q, self.fEx, self.fEy), interval=1, blit=False)

        draw()


if __name__ == '__main__':
    TE01 = pursue()
    TM01 = pursue()
    HE21 = pursue()
    TE01.te(0, 1)
    TM01.tm(0, 1)
    HE21.he(2, 1)

    mode1 = pursue()
    mode2 = pursue()

    mode1.te(0, 1)
    mode1 + HE21

    mode2.tm(0, 1)
    mode2 + HE21
    #
    # TE01.draw()
    # TM01.draw()
    # HE21.draw()
    #
    mode1.polarizer(45, 'off')
    mode1.movie()
    # mode2.draw()
    show()
