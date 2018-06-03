import numpy as np
import matplotlib.pyplot as plt
    
# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive
from pytriqs.gf import MeshBrillouinZone, Idx

from triqs_tprf.lattice_utils import k_space_path

# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_e_k_and_chi00_wk.h5'
    
    with HDFArchive(filename, 'r') as arch:
        e_k = arch['e_k']
        chi00_wk = arch['chi00_wk']

    e_k_interp = np.vectorize(lambda kx, ky, kz : e_k([kx, ky, kz])[0, 0].real)    
    
    G = np.array([0.0, 0.0, 0.0]) * 2.*np.pi
    X = np.array([0.5, 0.5, 0.0]) * 2.*np.pi
    M = np.array([0.5, 0.0, 0.0]) * 2.*np.pi
    
    paths = [(G, X), (X, M), (M, G)]

    k_vecs, k_plot, K_plot = k_space_path(paths)
    kx, ky, kz = k_vecs.T

    e_k_interp = e_k_interp(kx, ky, kz)

    plt.figure()
    
    plt.plot(k_plot, e_k_interp, '-')
    plt.axes().set_xticks(K_plot)
    plt.axes().set_xticklabels([r'$\Gamma$',r'$X$',r'$M$',r'$\Gamma$'])
    plt.ylabel(r'$\epsilon(\mathbf{k})$')
    plt.grid()

    plt.tight_layout()
    plt.savefig('figure_e_k_bandpath.pdf')

    # ------------------------------------------------------------------

    chi_k = chi00_wk[Idx(0), :][0, 0, 0, 0] # zero frequency, and scalar
    
    chi_k_interp = np.vectorize(lambda kx, ky, kz : chi_k([kx, ky, kz]).real)
    chi_k_interp = chi_k_interp(kx, ky, kz)
    
    plt.figure()
    
    plt.plot(k_plot, chi_k_interp, '-')
    plt.axes().set_xticks(K_plot)
    plt.axes().set_xticklabels([r'$\Gamma$',r'$X$',r'$M$',r'$\Gamma$'])
    plt.ylabel(r'$\chi_0(\mathbf{k}, i\omega_n = 0)$')
    plt.grid()

    plt.tight_layout()
    plt.savefig('figure_chi00_k_bandpath.pdf')    

    # ------------------------------------------------------------------
    plt.show()
    # ------------------------------------------------------------------
