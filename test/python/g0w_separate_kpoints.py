# ----------------------------------------------------------------------

import numpy as np

from triqs_tprf.lattice import g0w_sigma, g0w_dynamic_sigma

from triqs.gf import Gf, MeshReFreq, MeshBrZone
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)

def test_gw_separate_kpoints():
    """ Tests the various ways of giving a mesh to the g0w_sigma calculator against each other.
    Author: Yann in 't Veld (2023) """ 

    mu = 0.5
    g2 = 0.4
    wD = 0.1
    beta = 300.0
    nk = 10
    g2 = 0.1
    wD = 0.02

    wmin = -10.0
    wmax = 10.0
    nw = 2
    delta = 0.01

    # Construct kmesh with only Gamma point
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, np.diag(np.array([nk, nk, nk], dtype=int)))
    fmesh = MeshReFreq(wmin, wmax, nw)

    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        Enk.data[k.data_index,:] = 0.1*knorm**2.0

    print('--> bare interaction')
    # Some k-dependent non-physical interaction
    V_k = Gf(mesh=kmesh, target_shape=[1]*4)
    V_k.data[:] = 0.0
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        if(np.isclose(knorm, 0.0)): continue
        V_k.data[k.data_index,:,:] = knorm**2.0 + 1.0j / knorm

    print("--> dynamic interaction")
    W_fk = Gf(mesh=MeshProduct(fmesh,kmesh), target_shape=[1]*4)
    W_fk.data[:] = 0.0
    for f in fmesh:
        fii = f.data_index
        W_fk.data[fii,:] = ElectronPhononInteraction(f.value + 1.0j*delta, g2 ,wD) + V_k.data[:]

    print('--> g0w_sigma')
    print("full")
    sigma1_stat_k = g0w_sigma(mu, beta, Enk, V_k)
    sigma1_dyn_fk = g0w_dynamic_sigma(mu, beta, Enk, W_fk, V_k, delta)
    sigma1_fk = g0w_sigma(mu, beta, Enk, W_fk, V_k, delta)

    print("on a given kmesh")
    sigma2_stat_k = g0w_sigma(mu, beta, Enk, V_k, kmesh)
    sigma2_dyn_fk = g0w_dynamic_sigma(mu, beta, Enk, W_fk, V_k, delta, kmesh)
    sigma2_fk = g0w_sigma(mu, beta, Enk, W_fk, V_k, delta, kmesh)

    print("per k point")
    sigma3_stat_k = Gf(mesh=kmesh, target_shape=[1]*2)
    sigma3_dyn_fk = Gf(mesh=MeshProduct(fmesh,kmesh), target_shape=[1]*2)
    sigma3_fk = Gf(mesh=MeshProduct(fmesh,kmesh), target_shape=[1]*2)
    for k in kmesh:
        kii = k.data_index
        sigma3_stat_k.data[kii,:,:] = g0w_sigma(mu, beta, Enk, V_k, k.value)
        sigma3_dyn_fk.data[:,kii,:,:] = g0w_dynamic_sigma(mu, beta, Enk, W_fk, V_k, delta, k.value).data[:]
        sigma3_fk.data[:,kii,:,:] = g0w_sigma(mu, beta, Enk, W_fk, V_k, delta, k.value).data[:]
   
    print(sigma1_stat_k.data[0:10,0,0])

    print("--> compare static parts")
    np.testing.assert_array_almost_equal(sigma1_stat_k.data[:], sigma2_stat_k.data[:])
    np.testing.assert_array_almost_equal(sigma2_stat_k.data[:], sigma3_stat_k.data[:])
    np.testing.assert_array_almost_equal(sigma3_stat_k.data[:], sigma1_stat_k.data[:])

    print("--> compare dynamic parts")
    np.testing.assert_array_almost_equal(sigma1_dyn_fk.data[:], sigma2_dyn_fk.data[:])
    np.testing.assert_array_almost_equal(sigma2_dyn_fk.data[:], sigma3_dyn_fk.data[:])
    np.testing.assert_array_almost_equal(sigma3_dyn_fk.data[:], sigma1_dyn_fk.data[:])

    print("--> compare full sigma")
    np.testing.assert_array_almost_equal(sigma1_fk.data[:], sigma2_fk.data[:])
    np.testing.assert_array_almost_equal(sigma2_fk.data[:], sigma3_fk.data[:])
    np.testing.assert_array_almost_equal(sigma3_fk.data[:], sigma1_fk.data[:])

    print("--> compare full against sum")
    sigma1_ref_fk = sigma1_dyn_fk.copy()
    sigma2_ref_fk = sigma2_dyn_fk.copy()
    sigma3_ref_fk = sigma3_dyn_fk.copy()
   
    for f in fmesh:
        fii = f.data_index
        sigma1_ref_fk.data[fii,:] += sigma1_stat_k.data[:]
        sigma2_ref_fk.data[fii,:] += sigma2_stat_k.data[:]
        sigma3_ref_fk.data[fii,:] += sigma3_stat_k.data[:]

    np.testing.assert_array_almost_equal(sigma1_fk.data[:], sigma1_ref_fk.data[:])
    np.testing.assert_array_almost_equal(sigma2_fk.data[:], sigma2_ref_fk.data[:])
    np.testing.assert_array_almost_equal(sigma3_fk.data[:], sigma3_ref_fk.data[:])


if __name__ == "__main__":
    test_gw_separate_kpoints()
