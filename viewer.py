import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle
import numpy as np


class ModeShow:
    def __init__(self, mode_class):
        self.x = mode_class.xi
        self.a = mode_class.a
    
    def plot2d(self, fig, ax, mobj, **kwargs):
        self.plot_add_mode(ax, mobj, **kwargs)
        if kwargs.get('core') == 'on':
            self.plot_add_core(ax)
        if kwargs.get('vector') == 'on':
            self.plot_add_vector_field(ax, mobj) 
        if kwargs.get('cbar') == 'on':
            self.plot_add_colorbar(fig, ax)
        if kwargs.get('tickoff') == True:
            ax.set_xticks([])
            ax.set_yticks([])
        return fig, ax

    def plot_add_core(self, ax):
        core_edge = Circle((0,0), 
                        self.a, 
                        facecolor='none', 
                        edgecolor='w', 
                        linewidth=0.4, 
                        alpha=0.5, 
                        linestyle='-.')
        ax.add_patch(core_edge)
        return ax
    
    def plot_add_vector_field(self, ax, mobj):
        sp = 15
        Ex, Ey = mobj
        xv = self.x[::sp]
        Exv = Ex[::sp, ::sp]
        Eyv = Ey[::sp, ::sp]
        ax.quiver(xv, xv,
                Exv.real, Eyv.real,
                pivot="middle", 
                scale=0.3, 
                color='k', 
                linewidths=0.2, 
                headaxislength=3, 
                headlength=6, 
                headwidth=5)
        return ax
    
    def plot_add_colorbar(self, fig, ax, imnum=0):
        # get ScalarMappable object from ax
        im = ax.collections[imnum]
        cbar = fig.colorbar(im, format="%.1e")
        return cbar
    
    def plot_add_mode(self, ax, mobj, **kwargs):
        Ix, Iy = mobj.real**2 + mobj.imag**2 
        I = Ix + Iy # Intensity meshgrid
        if kwargs.get('type') == 'Ax':
            ax.pcolormesh(self.x, self.x, mobj[0].real, cmap=cm.twilight)
        elif kwargs.get('type') == 'Ay':
            ax.pcolormesh(self.x, self.x, mobj[1].real, cmap=cm.twilight)
        else:
            ax.pcolormesh(self.x, self.x, I, cmap=cm.nipy_spectral)
        return ax
