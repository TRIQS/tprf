# ----------------------------------------------------------------------

import numpy as np

from triqs_tprf.gw import g0w_sigma

from triqs.gf import Gf, MeshReFreq, MeshBrillouinZone
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

def test_gw_separate_kpoints():
    """ Tests the various ways of giving a mesh to the g0w_sigma calculator against each other.
    Author: Yann in 't Veld (2023) """ 

    mu = 0.5
    g2 = 0.4
    wD = 0.1
    beta = 300.0
    nk = 10

    # Construct kmesh with only Gamma point
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrillouinZone(bz, np.diag(np.array([nk, nk, nk], dtype=int)))

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

    print('--> g0w_sigma')
    print("full")
    sigma1_k = g0w_sigma(mu, beta, Enk, V_k)

    print("on a given kmesh")
    sigma2_k = g0w_sigma(mu, beta, Enk, V_k, kmesh)

    print("per k point")
    sigma3_k = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        kii = k.linear_index
        sigma3_k.data[kii,:,:] = g0w_sigma(mu, beta, Enk, V_k, k.value)
   
    print(sigma1_k.data[0:10,0,0])

    np.testing.assert_array_almost_equal(sigma1_k.data[:], sigma2_k.data[:])
    np.testing.assert_array_almost_equal(sigma2_k.data[:], sigma3_k.data[:])
    np.testing.assert_array_almost_equal(sigma3_k.data[:], sigma1_k.data[:])

if __name__ == "__main__":
    test_gw_separate_kpoints()
