# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.utility import mpi
from triqs.gf import Gf, MeshBrZone
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import gaussian_dos
from triqs_tprf.lattice import densdens_V_orb_to_band_basis
from triqs_tprf.lattice import densdens_V_pseudopotential

from triqs_tprf.lattice_utils import backfold_k

import scipy.integrate as integrate

# ----------------------------------------------------------------------

def get_analytical_pseudopotential_TF(a0, mstar, mu, eps):
    """Analytical expression for the Coulomb pseudo-potential for a Thomas-Fermi screened interaction.
    The derivation can be found in appendix A of the master thesis of Yann in 't Veld

    Note: There is a mistake in the expression in this source. The correct expression is:
    mu = N0 e^2 / (A epsilon k_F) I(2 pi e^2 / (A epsilon k_F) N0)"""

    A = a0*a0
    N0 = A*mstar / (2.0*np.pi)
    e2 = 14.399
    kF = np.sqrt(2.0 * mstar * mu)

    I = lambda x : integrate.quad(lambda theta:1.0/(np.sin(theta) + x), 0.0, np.pi)[0]
    mu_pseudo = N0 * e2 / (A * eps * kF) * I(2.0 * np.pi * e2 * N0 / (A * eps * kF))
    return mu_pseudo


def test_singleband():
    print("Single band")
    a0 = 2.0
    nk = 20
    sigma = 0.5
    mstar = 0.4 * 0.1314
    mu = 1.0

    UCA = a0*a0

    print("--> setup dispersion")
    bl = BravaisLattice(units=[(a0,0,0), (0,a0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, 1])
    eps_k = Gf(mesh=kmesh, target_shape=[1])
    for k in kmesh:
        kfolded = backfold_k(k.value, kmesh)
        knorm = np.linalg.norm(kfolded)
        eps_k[k] = knorm**2.0 / (2.0 * mstar)

    print('--> get dos')
    dos = gaussian_dos(eps_k, mu, sigma)

    print("--> setup Coulomb interaction")
    eps = 1.0
    e2 = 14.399

    W_q = Gf(mesh=kmesh, target_shape=[1]*4)
    W_q.data[:] = 0.0
    for k in kmesh:
        qfolded = backfold_k(k.value, kmesh)
        qnorm = np.linalg.norm(qfolded)

        invV = UCA * qnorm * eps / (2.0 * np.pi * e2)
        W_q[k] = 1.0 / (invV + 2.0*dos)

    psi_k = Gf(mesh=kmesh, target_shape=[1]*2)
    psi_k.data[:] = 1.0

    print("--> densdens_V_pseudopotential")
    Vb_kkp = densdens_V_orb_to_band_basis(W_q, psi_k)
    mu_pseudo = densdens_V_pseudopotential(eps_k, mu, sigma, Vb_kkp)
    print("num:", mu_pseudo.real)

    mu_pseudo_ref = get_analytical_pseudopotential_TF(a0, mstar, mu, eps)
    print("ana:", mu_pseudo_ref)

    assert np.isclose(mu_pseudo, mu_pseudo_ref, rtol=1e-2)


if __name__ == "__main__":
    test_singleband()
