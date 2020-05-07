
import numpy as np
import matplotlib.pyplot as plt
    
# ----------------------------------------------------------------------

from scipy.interpolate import griddata

# ----------------------------------------------------------------------

from h5 import HDFArchive
from triqs.gf import MeshBrillouinZone, Idx

# ----------------------------------------------------------------------
def plot_chi(chi, title):

    chi_k = chi[Idx(0), :][0, 0, 0, 0] # zero frequency, and scalar
    
    k = np.linspace(-0.75, 0.75, num=100) * 2. * np.pi
    Kx, Ky = np.meshgrid(k, k)

    interp = np.vectorize(lambda kx, ky : chi_k([kx, ky, 0]).real)    
    interp = interp(Kx, Ky)

    plt.imshow(
        interp, cmap=plt.get_cmap('magma'),
        extent=(k.min(), k.max(), k.min(), k.max()),
        origin='lower', vmin=0, vmax=interp.max())

    plt.title(title)
    plt.xlabel(r'$\frac{a}{2\pi} \cdot k_x$')
    plt.ylabel(r'$\frac{a}{2\pi} \cdot k_y$')
    plt.colorbar()

# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_e_k_and_chi00_wk.h5'
    
    with HDFArchive(filename, 'r') as arch:
        chi00_wk = arch['chi00_wk']
        chi00_wk_analytic = arch['chi00_wk_analytic']

    plt.figure(figsize=(8, 3))

    subp = [1, 2, 1]

    plt.subplot(*subp); subp[-1] += 1
    plot_chi(chi00_wk,
             title=r'$\chi^{00}_{S_z S_z}(\mathbf{q}, \omega=0)$ (imtime)')

    plt.subplot(*subp); subp[-1] += 1
    plot_chi(chi00_wk_analytic,
             title=r'$\chi^{00}_{S_z S_z}(\mathbf{q}, \omega=0)$ (analytic)')
    
    plt.tight_layout()
    plt.savefig('figure_lindhardt_SzSz_square_latt.pdf')
    plt.show()

