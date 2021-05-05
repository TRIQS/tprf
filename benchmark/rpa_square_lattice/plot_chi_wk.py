
import numpy as np
import matplotlib.pyplot as plt
    
# ----------------------------------------------------------------------

from scipy.interpolate import griddata

# ----------------------------------------------------------------------

from h5 import HDFArchive
from triqs.gf import MeshBrZone, Idx
from triqs.lattice import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

from triqs_tprf.lattice_utils import k_space_path
from triqs_tprf.lattice_utils import get_relative_k_from_absolute

# ----------------------------------------------------------------------
def get_chi_SzSz(chi):

    # -- Contract Sz Sz response
    Sx = 0.5 * np.array([[0., 1.], [1., 0]])
    Sy = 0.5 * np.array([[0., -1.j], [1.j, 0]])
    Sz = np.diag([+0.5, -0.5])
    N = np.eye(2)

    #Op1, Op2 = Sx, Sx
    #Op1, Op2 = Sy, Sy
    Op1, Op2 = Sz, Sz
    #Op1, Op2 = N, N
    
    chi_SzSz = chi[0, 0, 0, 0].copy()
    chi_SzSz.data[:] = np.einsum('wqabcd,ab,cd->wq', chi.data, Op1, Op2)[:, :]
    chi_SzSz = chi_SzSz[Idx(0), :]
    return chi_SzSz

# ----------------------------------------------------------------------
def plot_chi(chi, label=None):

    chi_SzSz = get_chi_SzSz(chi)

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

    chi_SzSz = get_chi_SzSz(chi)

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

    # -- band path plot
    
    plt.figure(figsize=(6, 5))

    plot_chi_1D(chi00_wk)

    for U, chi_wk in zip(U_vec, chi_wk_vec):
        plot_chi_1D(chi_wk, label='U=%2.2f' % U)

    plt.legend(loc='best')
    plt.savefig('figure_RPA_SzSz_square_latt_bandpath.pdf')

    # -- Extrapolate critical U_c
    
    plt.figure(figsize=(6, 5))

    inv_chi_vec = []
    
    for U, chi_wk in zip(U_vec, chi_wk_vec):
        chi_SzSz = get_chi_SzSz(chi_wk)
        inv_chi = np.min(1./chi_SzSz.data.real)
        plt.plot(U, inv_chi, 'kx')
        inv_chi_vec.append(inv_chi)

    p = np.polyfit(U_vec, inv_chi_vec, 1)
    Uc = np.roots(p)
    plt.plot(Uc, 0, 'rx', label=r'$U_c = %2.2f$' % Uc)
    plt.legend()
    plt.xlabel(r'$U$')
    plt.ylabel(r'min[$\chi_{RPA}^{-1}$]')
    
    # -- k-xy plane chi plot

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

