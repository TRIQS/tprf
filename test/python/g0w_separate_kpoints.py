# ----------------------------------------------------------------------

import numpy as np

#from triqs_tprf.gw import g0w_sigma, g0w_dyn_sigma
from triqs_tprf.lattice import g0w_sigma, g0w_dyn_sigma

from triqs.gf import Gf, MeshReFreq, MeshBrillouinZone
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
    kmesh = MeshBrillouinZone(bz, np.diag(np.array([nk, nk, nk], dtype=int)))
    fmesh = MeshReFreq(wmin, wmax, nw)

    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        Enk.data[k.linear_index,:] = 0.1*knorm**2.0

    print('--> bare interaction')
    # Some k-dependent non-physical interaction
    V_k = Gf(mesh=kmesh, target_shape=[1]*4)
    V_k.data[:] = 0.0
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        if(np.isclose(knorm, 0.0)): continue
        V_k.data[k.linear_index,:,:] = knorm**2.0 + 1.0j / knorm

    print("--> dynamic interaction")
    W_fk = Gf(mesh=MeshProduct(fmesh,kmesh), target_shape=[1]*4)
    W_fk.data[:] = 0.0
    for f in fmesh:
        fii = f.linear_index
        W_fk.data[fii,:] = ElectronPhononInteraction(f.value + 1.0j*delta, g2 ,wD) + V_k.data[:]

    print('--> g0w_sigma')
    print("full")
    sigma1_k = g0w_sigma(mu, beta, Enk, V_k)
    sigma1_fk = g0w_dyn_sigma(mu, beta, Enk, W_fk, V_k, delta)

    print("on a given kmesh")
    sigma2_k = g0w_sigma(mu, beta, Enk, V_k, kmesh)
    sigma2_fk = g0w_dyn_sigma(mu, beta, Enk, W_fk, V_k, delta, kmesh)

    print("per k point")
    sigma3_k = Gf(mesh=kmesh, target_shape=[1]*2)
    sigma3_fk = Gf(mesh=MeshProduct(fmesh,kmesh), target_shape=[1]*2)
    for k in kmesh:
        kii = k.linear_index
        sigma3_k.data[kii,:,:] = g0w_sigma(mu, beta, Enk, V_k, k.value)
        sigma3_fk.data[:,kii,:,:] = g0w_dyn_sigma(mu, beta, Enk, W_fk, V_k, delta, k.value).data[:]
   
    print(sigma1_k.data[0:10,0,0])

    print("--> compare static parts")
    np.testing.assert_array_almost_equal(sigma1_k.data[:], sigma2_k.data[:])
    np.testing.assert_array_almost_equal(sigma2_k.data[:], sigma3_k.data[:])
    np.testing.assert_array_almost_equal(sigma3_k.data[:], sigma1_k.data[:])

    print("--> compare dynamic parts")
    np.testing.assert_array_almost_equal(sigma1_fk.data[:], sigma2_fk.data[:])
    np.testing.assert_array_almost_equal(sigma2_fk.data[:], sigma3_fk.data[:])
    np.testing.assert_array_almost_equal(sigma3_fk.data[:], sigma1_fk.data[:])

if __name__ == "__main__":
    test_gw_separate_kpoints()
