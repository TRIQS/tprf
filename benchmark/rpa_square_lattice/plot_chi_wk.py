
import numpy as np
import matplotlib.pyplot as plt
    
# ----------------------------------------------------------------------

from scipy.interpolate import griddata

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive
from pytriqs.gf import MeshBrillouinZone, Idx
from pytriqs.lattice import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

from triqs_tprf.lattice_utils import k_space_path
from triqs_tprf.lattice_utils import get_relative_k_from_absolute

# ----------------------------------------------------------------------
def plot_chi(chi, label=None):

    # -- Contract Sz Sz response
    Sz = np.diag([+0.5, -0.5])
    chi_SzSz = chi[0, 0, 0, 0].copy()
    chi_SzSz.data[:] = np.einsum('wqabcd,ab,cd->wq', chi.data, Sz, Sz)[:, :]
    chi_SzSz = chi_SzSz[Idx(0), :]

    k = np.linspace(-0.75, 0.75, num=100) * 2.*np.pi
    Kx, Ky = np.meshgrid(k, k)
    k_vecs = np.vstack((Kx.flatten(), Ky.flatten(), 0*Kx.flatten())).T

    chi_interp = np.zeros(
        [k_vecs.shape[0]] + list(chi_SzSz.target_shape), dtype=np.complex)

    for kidx, (kx, ky, kz) in enumerate(k_vecs):
        chi_interp[kidx] = chi_SzSz((kx, ky, kz))

    chi_interp = chi_interp.real.reshape(Kx.shape)

    plt.imshow(
        chi_interp, cmap=plt.get_cmap('magma'),
        extent=(k.min(), k.max(), k.min(), k.max()),
        origin='lower',
        vmin=0, vmax=chi_interp.max(),
        )

    plt.title(label)
    #plt.title(r'Lindhardt spin-response $\chi^{00}_{S_z S_z}(\mathbf{q}, \omega=0)$')
    plt.xlabel(r'$\frac{a}{2\pi} \cdot k_x$')
    plt.ylabel(r'$\frac{a}{2\pi} \cdot k_y$')
    plt.colorbar()

# ----------------------------------------------------------------------
def plot_chi_1D(chi, label=None):

    # -- Contract Sz Sz response
    Sz = np.diag([+0.5, -0.5])
    chi_SzSz = chi[0, 0, 0, 0].copy()
    chi_SzSz.data[:] = np.einsum('wqabcd,ab,cd->wq', chi.data, Sz, Sz)[:, :]
    chi_SzSz = chi_SzSz[Idx(0), :]

    bz = chi_SzSz.mesh.domain

    G = np.array([0., 0., 0.]) * 2.*np.pi
    X = np.array([0.5, 0., 0.]) * 2.*np.pi
    M = np.array([0.5, 0.5, 0.]) * 2.*np.pi

    paths = [ (G, M), (M, X), (X, G) ]
    k_vecs, k_plot, K_plot = k_space_path(paths)
    k_vecs = get_relative_k_from_absolute(k_vecs, bz.units())
    k_vecs *= 2.*np.pi

    chi_interp = np.zeros(
        [k_vecs.shape[0]] + list(chi_SzSz.target_shape), dtype=np.complex)

    for kidx, (kx, ky, kz) in enumerate(k_vecs):
        chi_interp[kidx] = chi_SzSz((kx, ky, kz))

    chi_interp = chi_interp.real

    plt.plot(k_plot, chi_interp, label=label)
    
    plt.grid()
    plt.axes().set_xticks(K_plot)
    plt.xlim([K_plot.min(), K_plot.max()])
    plt.axes().set_xticklabels([r'$\Gamma$',r'$M$',r'$X$',r'$\Gamma$'])
    
    plt.title(r'RPA spin-response $\chi_{S_z S_z}(\mathbf{q}, \omega=0)$')
    
# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_ek_and_chi_wk.h5'
    
    with HDFArchive(filename, 'r') as arch:
        chi00_wk = arch['chi00_wk']
        chi_wk_vec = arch['chi_wk_vec']
        U_vec = arch['U_vec']

    plt.figure(figsize=(6, 5))

    plot_chi_1D(chi00_wk)

    for U, chi_wk in zip(U_vec, chi_wk_vec):
        plot_chi_1D(chi_wk, label='U=%2.2f' % U)

    plt.legend(loc='best')
    plt.savefig('figure_RPA_SzSz_square_latt_bandpath.pdf')

    plt.figure(figsize=(8, 3))

    subp = [1, 2, 1]

    plt.subplot(*subp); subp[-1] += 1
    plot_chi(
        chi00_wk,
        label=r'Lindhardt spin-response $\chi^{00}_{S_z S_z}(\mathbf{q}, \omega=0)$')

    plt.subplot(*subp); subp[-1] += 1
    plot_chi(
        chi_wk_vec[-1],
        label=r'RPA U=%2.2f $\chi^{00}_{S_z S_z}(\mathbf{q}, \omega=0)$' % U_vec[-1])
    
    plt.tight_layout()
    plt.savefig('figure_RPA_SzSz_square_latt_plane.pdf')
    plt.show()

