# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk, lattice_dyson_g_w

from triqs_tprf.gw import bubble_PI_wk
from triqs_tprf.gw import dynamical_screened_interaction_W_wk
from triqs_tprf.gw import gw_sigma_wk

from triqs.gf import Gf, MeshImFreq, MeshReFreq
from triqs.gf.mesh_product import MeshProduct

# ----------------------------------------------------------------------

def test_gf_Matsubara():
    nw = 100
    nk = 8
    norb = 2
    beta = 10.0
    V = 5.0
    mu = 0.0
    
    # Tight-binding
    t = -1.0 * np.eye(norb)
    
    t_r = TBLattice(
        units = [(1, 0, 0)],
        hopping = {
            (+1,) : t,
            (-1,) : t,
            },
        orbital_positions = [(0,0,0)]*norb,
        )
    
    kmesh = t_r.get_kmesh(n_k=(nk, 1, 1))
    e_k = t_r.fourier(kmesh)
    
    kmesh = e_k.mesh
    wmesh = MeshImFreq(beta, 'Fermion', nw)

    # Test setup bare Gf
    print("  -> bare Matsubara Gf")
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

    g0_wk_ref = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=[norb]*2)
    for w in wmesh:
        for k in kmesh:
            g0_wk_ref[w,k] = np.linalg.inv( (w.value + mu)*np.eye(norb) - e_k[k] )

    np.testing.assert_array_almost_equal(g0_wk.data[:], g0_wk_ref.data[:])

    # Get self-energy
    V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_k.data[:] = V
    PI_wk = bubble_PI_wk(g0_wk)
    W_wk = dynamical_screened_interaction_W_wk(PI_wk, V_k)
    sigma_wk = gw_sigma_wk(W_wk, g0_wk)
    
    # Construct a kmesh-independent self-energy
    sigma_w = Gf(mesh=wmesh, target_shape=[norb]*2)
    sigma_w.data[:] = sigma_wk.data[:,0]


    # lattice_dyson_g_wk, input sigma_wk
    print("  -> lattice_dyson_g_wk, sigma_wk")
    g_wk = lattice_dyson_g_wk(mu=mu, e_k=e_k, sigma_wk=sigma_wk)
    g_wk_ref = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=[norb]*2)
    for w in wmesh:
        for k in kmesh:
            g_wk_ref[w,k] = np.linalg.inv( (w.value + mu)*np.eye(norb) - e_k[k] - sigma_wk[w,k] )
    
    np.testing.assert_array_almost_equal(g_wk.data[:], g_wk_ref.data[:])

    
    # lattice_dyson_g_wk, input sigma_w
    print("  -> lattice_dyson_g_wk, sigma_w")
    g_wk_2 = lattice_dyson_g_wk(mu=mu, e_k=e_k, sigma_w=sigma_w)
    g_wk_2_ref = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=[norb]*2)
    for w in wmesh:
        for k in kmesh:
            g_wk_2_ref[w,k] = np.linalg.inv( (w.value + mu)*np.eye(norb) - e_k[k] - sigma_w[w] )
    
    np.testing.assert_array_almost_equal(g_wk_2.data[:], g_wk_2_ref.data[:])


    # lattice_dyson_g_w, input sigma_w
    print("  -> lattice_dyson_g_w, sigma_w")
    g_w = lattice_dyson_g_w(mu=mu, e_k=e_k, sigma_w=sigma_w)
    
    g_w_ref = Gf(mesh=wmesh, target_shape=[norb]*2)
    g_w_ref.data[:] = 0.0
    for w in wmesh:
        wi = w.linear_index
        for k in kmesh:
            g_w_ref.data[wi,:] += np.linalg.inv( (w.value + mu)*np.eye(norb) - e_k[k] - sigma_w[w] )
    g_w_ref.data[:] /= len(kmesh)

    #g_w_ref = Gf(mesh=wmesh, target_shape=[norb]*2)
    #g_w_ref.data[:] = 0.0
    #for w in wmesh:
    #    wi = w.linear_index
    #    g_w_ref.data[wi,:] += np.sum(g_wk_2_ref.data[wi,:], axis=0)
    #g_w_ref.data[:] /= len(kmesh)

    np.testing.assert_array_almost_equal(g_w.data[:], g_w_ref.data[:])


if __name__ == "__main__":
    test_gf_Matsubara()




