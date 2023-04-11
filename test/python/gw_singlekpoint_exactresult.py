# ----------------------------------------------------------------------

import numpy as np

from triqs_tprf.lattice import lattice_dyson_g0_wk
#from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.gw import gw_sigma, g0w_sigma

from triqs.gf import Gf, MeshImFreq, MeshReFreq, MeshBrZone
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

def nF(ww, beta):
    """Fermi-Dirac distribution function"""
    return 0.5 - 0.5 * np.tanh(0.5 * beta * ww)

def nB(ww, beta):
    """Bose-Einstein distribution function"""
    return 1.0 / np.expm1( beta * ww )

def ExactSigma0D(iw, beta, g2, wD, E):
    """Exact result for the 0-dimensional GW self-energy for an electron-phonon interaction with dispersionless
    dispersion and scalar electron-phonon coupling, using a single k-point at gamma."""
    fact1 = (nB(wD, beta) + nF(E, beta)) / (iw + wD - E)
    fact2 = (nB(wD, beta) + 1.0 - nF(E, beta)) / (iw - wD - E)
    return g2 * (fact1 + fact2)

def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)

def test_gw_sigma_against_exact_Matsubara():
    """ Tests the Matsubara frequency axis GW implementation by comparing to an exact analytic result.
    This result was found for a calculation on a single k-point, with a simple electron-phonon propagator.
    Author: Yann in 't Veld (2023) """ 
    
    print("== Matsubara axis ==")

    mu = 0.5
    g2 = 0.4
    wD = 0.1
    beta = 300.0
    nw = 500

    # Construct kmesh with only Gamma point
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, np.diag(np.array([1, 1,1], dtype=int)))
    wmesh = MeshImFreq(beta, 'Fermion', nw)
    numesh = MeshImFreq(beta, 'Boson', nw)

    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    Enk.data[:] = 0.0
    g0_fk = lattice_dyson_g0_wk(mu, Enk, wmesh)

    print('--> bare electron-phonon interaction')
    I_phon_wk = Gf(mesh=MeshProduct(numesh, kmesh), target_shape=[1]*4)
    for nu in numesh:
        nuii = nu.data_index
        I_phon_wk.data[nuii,:] = ElectronPhononInteraction(nu.value, g2 ,wD)

    print('--> gw_sigma')
    sigma_wk = gw_sigma(I_phon_wk, g0_fk)

    sigma_ref_wk = Gf(mesh=sigma_wk.mesh, target_shape=sigma_wk.target_shape)
    for f in wmesh:
        fii = f.data_index
        sigma_ref_wk.data[fii,:] = ExactSigma0D(f.value, beta, g2, wD, Enk.data[0,0,0]-mu)

    np.testing.assert_array_almost_equal(sigma_wk.data[:], sigma_ref_wk.data[:], decimal=1e-6)


def test_gw_sigma_against_exact_realfreq():
    """ Tests thereal-frequency axis GW implementation by comparing to an exact analytic result.
    This result was found for a calculation on a single k-point, with a simple electron-phonon propagator.
    Author: Yann in 't Veld (2023) """ 

    print("== Real-freq axis ==")

    mu = 0.5
    g2 = 0.4
    wD = 0.1
    beta = 300.0
    nw = 5000
    wmin = -2.0
    wmax = 2.0
    eta = 1e-3

    # Construct kmesh with only Gamma point
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, np.diag(np.array([1, 1,1], dtype=int)))
    fmesh = MeshReFreq(wmin, wmax, nw)

    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    Enk.data[:] = 0.0

    print('--> bare electron-phonon interaction')
    I_phon_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=[1]*4)
    for f in fmesh:
        fii = f.data_index
        I_phon_fk.data[fii,:] = ElectronPhononInteraction(f.value + 1.0j*eta, g2 ,wD)

    I_phon_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_phon_k.data[:] = 0.0

    print('--> gw_sigma')
    sigma_fk = g0w_sigma(mu, beta, Enk, I_phon_fk, I_phon_k, eta)

    sigma_ref_fk = Gf(mesh=sigma_fk.mesh, target_shape=sigma_fk.target_shape)
    for f in fmesh:
        fii = f.data_index
        sigma_ref_fk.data[fii,:] = ExactSigma0D(f.value + 1.0j*eta, beta, g2, wD, Enk.data[0,0,0]-mu)

    # Rather large tolerance needed. Probably due to the numerical spectral function of I not being a
    # true delta-function, and thus deviating from the analytic result.
    assert np.allclose(sigma_ref_fk.data[:], sigma_fk.data[:], rtol=1.0)


if __name__ == "__main__":
    test_gw_sigma_against_exact_Matsubara()
    test_gw_sigma_against_exact_realfreq()
