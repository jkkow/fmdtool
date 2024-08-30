import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle


class ModeShow:
    def __init__(self, mode_class):
        self.x = mode_class.xi
        self.a = mode_class.a

    def plot2d(self, fig, axe, mobj, **kwargs):
        self.add_intensity_to_axe(axe, mobj, **kwargs)
        if kwargs.get("core") == "on":
            self.add_core_to_axe(axe)
        if kwargs.get("vector") == "on":
            self.add_vector_field_to_axe(axe, mobj)
        if kwargs.get("cbar") == "on":
            self.add_colorbar_to_axe(fig, axe)
        if kwargs.get("tickoff") == True:
            axe.set_xticks([])
            axe.set_yticks([])
        return fig, axe

    def add_intensity_to_axe(self, axe, mobj, **kwargs):
        Ix, Iy = mobj.real**2 + mobj.imag**2
        I = Ix + Iy  # Intensity meshgrid
        if kwargs.get("type") == "Ax":
            axe.pcolormesh(self.x, self.x, mobj[0].real, cmap=cm.twilight)
        elif kwargs.get("type") == "Ay":
            axe.pcolormesh(self.x, self.x, mobj[1].real, cmap=cm.twilight)
        else:
            axe.pcolormesh(self.x, self.x, I, cmap=cm.nipy_spectral)
        return axe

    def add_core_to_axe(self, axe):
        core_edge = Circle(
            (0, 0),
            self.a,
            facecolor="none",
            edgecolor="w",
            linewidth=0.4,
            alpha=0.5,
            linestyle="-.",
        )
        axe.add_patch(core_edge)
        return axe

    def add_vector_field_to_axe(self, axe, mobj):
        sp = 15
        Ex, Ey = mobj
        xv = self.x[::sp]
        Exv = Ex[::sp, ::sp]
        Eyv = Ey[::sp, ::sp]
        axe.quiver(
            xv,
            xv,
            Exv.real,
            Eyv.real,
            pivot="middle",
            scale=0.3,
            color="k",
            linewidths=0.2,
            headaxislength=3,
            headlength=6,
            headwidth=5,
        )
        return axe

    def add_colorbar_to_axe(self, fig, axe, imnum=0):
        # get ScalarMappable object from ax
        im = axe.collections[imnum]
        cbar = fig.colorbar(im, format="%.1e")
        return cbar
