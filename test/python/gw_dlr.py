# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.lattice import dlr_on_imfreq, dlr_on_imtime

from triqs_tprf.gw import bubble_PI_wk
from triqs_tprf.gw import dynamical_screened_interaction_W
from triqs_tprf.gw import gw_sigma, fock_sigma

from triqs.gf import Gf, Idx
from triqs.gf import MeshImFreq, MeshImTime, MeshDLRImFreq, MeshDLRImTime
from triqs.gf.mesh_product import MeshProduct

from triqs.gf.gf_factories import make_gf_dlr

# ----------------------------------------------------------------------

def dlr_wk_on_imfreq_wk(g_Dwk, wmesh, g_k=None):
    DLRwmesh, kmesh = g_Dwk.mesh.components
    g_wk =Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=g_Dwk.target_shape)
    g_wk.data[:] = 0.0

    if(g_k is None):
        g_k = Gf(mesh=kmesh, target_shape=g_Dwk.target_shape)
        g_k.data[:] = 0.0

    for k in kmesh:
        g_Dw = Gf(mesh=DLRwmesh, target_shape=g_Dwk.target_shape)
        g_Dw.data[:] = g_Dwk.data[:,k.data_index,:] - g_k.data[k.data_index,:]
        g_Dc = make_gf_dlr(g_Dw)
        g_wk[:,k] = dlr_on_imfreq(g_Dc, wmesh)
        g_wk.data[:,k.data_index,:] += g_k.data[k.data_index,:]

    return g_wk

def compare_g_Dwk_and_g_wk(g_Dwk, g_wk, g_k=None, decimal=7):
    wmesh, kmesh = g_wk.mesh.components
    DLRwmesh = g_Dwk.mesh.components[0]

    g_ref_wk = dlr_wk_on_imfreq_wk(g_Dwk, wmesh, g_k)

    np.testing.assert_array_almost_equal(g_wk.data[:], g_ref_wk.data[:], decimal=decimal)


def test_gw_sigma_functions():
    nk = 8
    norb = 1
    beta = 10.0
    V = 1.0
    mu = 0.1
    dlreta = 1e-10
    dlrwcut = 100.0
    
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

    wmesh_dlr = MeshDLRImFreq(beta, 'Fermion', dlrwcut, dlreta)

    nw = int(dlrwcut * beta / (2.0 * np.pi) - 0.5)
    wmesh = MeshImFreq(beta, 'Fermion', nw)

    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
    g0_Dwk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh_dlr)
    
    V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_k.data[:] = V
    
    print('--> pi_bubble')
    PI_wk = bubble_PI_wk(g0_wk)
    PI_Dwk = bubble_PI_wk(g0_Dwk)

    compare_g_Dwk_and_g_wk(PI_Dwk, PI_wk)
   
    print('--> gw_sigma (static)')
    sigma_k = gw_sigma(V_k, g0_wk)
    sigma_dlr_k = gw_sigma(V_k, g0_Dwk)

    np.testing.assert_array_almost_equal(sigma_k.data[:], sigma_dlr_k.data[:])

    print('--> dynamical_screened_interaction_W')
    Wr_full_wk = dynamical_screened_interaction_W(PI_wk, V_k)
    Wr_full_Dwk = dynamical_screened_interaction_W(PI_Dwk, V_k)

    compare_g_Dwk_and_g_wk(Wr_full_Dwk, Wr_full_wk, V_k, decimal=7)

    print('--> gw_sigma (dynamic)')
    sigma_wk = gw_sigma(Wr_full_wk, g0_wk)
    sigma_Dwk = gw_sigma(Wr_full_Dwk, V_k, g0_Dwk)

    compare_g_Dwk_and_g_wk(sigma_Dwk, sigma_wk, sigma_k, decimal=6)

    print('--> lattice_dyson_g_wk')
    g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
    g_Dwk = lattice_dyson_g_wk(mu, e_k, sigma_Dwk)

    compare_g_Dwk_and_g_wk(g_Dwk, g_wk)


if __name__ == "__main__":
    test_gw_sigma_functions()

