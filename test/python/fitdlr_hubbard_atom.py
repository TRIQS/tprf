
import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import Gf, MeshImTime, MeshImFreq
from triqs.operators import c, c_dag, n

from triqs.gf import inverse, iOmega_n, Fourier
from triqs.gf.gf_factories import make_gf_imtime, make_gf_imfreq, make_gf_dlr_imtime

from triqs.atom_diag import AtomDiag, atomic_g_iw, atomic_g_tau

# ----------------------------------------------------------------------

from triqs_tprf.fitdlr import fitdlr
from triqs_tprf.fitdlr import BlockSymmetrizer

# ----------------------------------------------------------------------

def get_H(U=1., mu=0.0, B=0.05):
    H = U * (n('0', 0) - 0.5)*(n('0', 1) - 0.5) - mu * (n('0', 0) + n('0', 1)) + B * (n('0', 0) - n('0', 1))
    return H


def get_gfs(H, beta, fundamental_operators, niw=128):

    ntau = 6 * niw + 1

    gf_struct = [('0', 2)]
    fops = [ list(o)[0][0][-1][-1] for o in fundamental_operators ]

    ad = AtomDiag(H, fops)
    G_iw = atomic_g_iw(ad, beta, gf_struct, niw)['0']
    G_tau = atomic_g_tau(ad, beta, gf_struct, ntau)['0']
    
    return G_tau, G_iw
    
        
def test_fit(verbose=False):

    np.random.seed(seed=1337)
    
    xi = -1.
    beta = 20.0

    H0 = get_H(U=0.)
    H = get_H()
    H_int = H - H0

    fundamental_operators = [c('0', 0), c('0', 1)]

    G0_tau, G0_iw = get_gfs(H0, beta, fundamental_operators)
    G_tau, G_iw = get_gfs(H, beta, fundamental_operators)

    Sigma_iw = G_iw.copy()
    Sigma_iw << inverse(G0_iw) - inverse(G_iw)

    # -- Fit greens function with DLR expansion

    G_tau_org = G_tau.copy()
    tol = 1e-5
    G_tau.data[:] += np.random.normal(scale=tol, size=G_tau.data.shape)

    from triqs.gf.meshes import MeshDLR
    cmesh = MeshDLR(beta, 'Fermion', w_max=1., eps=1e-6)
        
    block_mat = np.array([[1, 0], [0, 2]])
    sym = BlockSymmetrizer(len(cmesh), block_mat)

    opt = dict(
        discontinuity=True,
        density=True,
        realvalued=True,
        ftol=1e-6,
        symmetrizer=sym,
        verbose=verbose,
        )
    
    G_c, sol = fitdlr(cmesh, G_tau, H, fundamental_operators, **opt)
    G_xaa_sym = G_c.data
    
    G_l_sym = make_gf_dlr_imtime(G_c)
    G_laa_sym = G_l_sym.data

    G_l = G_l_sym
    G_xaa = G_xaa_sym
    G_laa = G_laa_sym

    tau_i = np.array([float(t) for t in G_tau.mesh])
    
    G_tau_fit = make_gf_imtime(G_c, len(tau_i))
    G_tau_diff = np.max(np.abs(G_tau_org.data - G_tau_fit.data))
    
    iwn = np.array([complex(w) for w in G_iw.mesh])
    G_iw_fit = make_gf_imfreq(G_c, len(G_iw.mesh)//2)
    G_iw_diff = np.max(np.abs(G_iw.data - G_iw_fit.data))
    
    Sigma_iw_fit = G_iw_fit.copy()
    Sigma_iw_fit << inverse(G0_iw) - inverse(G_iw_fit)
    Sigma_iw_diff = np.max(np.abs(Sigma_iw.data - Sigma_iw_fit.data))

    if verbose:
        print(f'G_tau_diff = {G_tau_diff:1.1E}')
        print(f'G_iw_diff = {G_iw_diff:1.1E}')
        print(f'Sigma_iw_diff = {Sigma_iw_diff:1.1E}')

    assert( G_tau_diff < 1e-4 )
    assert( G_iw_diff < 1e-4 )
    assert( Sigma_iw_diff < 1e-3 )

    if verbose:
        tau_l = np.array([ float(t) for t in G_l.mesh ])

        from triqs.plot.mpl_interface import oplot, oplotr, oploti, plt

        plt.figure(figsize=(12, 12))
        subp = [4, 2, 1]

        plt.subplot(*subp); subp[-1] += 1
        oplotr(G0_tau)
        oploti(G0_tau, '--')
        plt.ylabel(r'$G_0(\tau)$')

        plt.subplot(*subp); subp[-1] += 1
        oplotr(-G_tau, alpha=0.25)
        oploti(-G_tau, '--', alpha=0.25)
        oplotr(-G_tau_fit)
        oploti(-G_tau_fit, '--')
        for i, j in itertools.product(range(2), repeat=2):
            plt.plot(tau_l, -G_laa[:, i, j].real, '+')
            plt.plot(tau_l, -G_laa[:, i, j].imag, 'x')
        plt.ylabel(r'$G(\tau)$')
        plt.semilogy([], [])

        plt.legend(loc='best')

        plt.subplot(*subp); subp[-1] += 1
        oplot(G0_iw, label='G0')
        plt.subplot(*subp); subp[-1] += 1
        oplot(G_iw, label='G')

        plt.subplot(*subp); subp[-1] += 1
        oplotr(Sigma_iw_fit, '.')
        oplotr(Sigma_iw)
        plt.ylabel(r'Re[$\Sigma$]')

        plt.subplot(*subp); subp[-1] += 1
        oploti(Sigma_iw_fit, '.')
        oploti(Sigma_iw)
        plt.ylabel(r'Im[$\Sigma$]')

        plt.subplot(*subp); subp[-1] += 1
        oplotr(Sigma_iw - Sigma_iw_fit, '.')
        plt.ylabel(r'Re[$\Delta \Sigma$]')

        plt.subplot(*subp); subp[-1] += 1
        oploti(Sigma_iw - Sigma_iw_fit, '.')
        plt.ylabel(r'Im[$\Delta \Sigma$]')
        
        plt.tight_layout()
        plt.show()
            
    
if __name__ == '__main__':
    test_fit(verbose=False)
