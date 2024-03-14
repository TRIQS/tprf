# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k, add_dynamic_and_static
from triqs_tprf.lattice import g_wk_to_g_mwk

from triqs_tprf.gw import bubble_PI_wk
from triqs_tprf.gw import dynamical_screened_interaction_W

from triqs.gf import Gf, MeshImFreq, Idx, MeshImTime, MeshBrZone
from triqs.gf.meshes import MeshDLRImFreq, MeshReFreq
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice
from triqs.gf.gf_factories import make_gf_dlr

# ----------------------------------------------------------------------

def test_g_wk_to_g_mwk():
    print("== DLR G(w,k) to G(-w,k) ==")
    beta = 10.0
    nk = 5

    lamb = 10.
    eps = 1e-8

    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, nk])
    wmesh = MeshDLRImFreq(beta, 'Fermion', lamb, eps)

    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        Enk.data[k.data_index,:] = 0.1 * knorm**2.0

    print("--> g_wk_to_g_mwk")
    g_wk = lattice_dyson_g0_wk(0.2, Enk, wmesh)
    g_mwk = g_wk_to_g_mwk(g_wk)

    print("--> check")
    for k in kmesh:
        g_c = make_gf_dlr(g_wk[:,k])
        g_mc = make_gf_dlr(g_mwk[:,k])

        for w in wmesh:
            assert np.allclose(g_c(w), g_mc(-w))


def test_split_into_dynamic_and_constant():
    print("== Split in to dynamic and constant ==")
    nw = 1000
    nk = 8
    norb = 1
    beta = 10.0
    V = 1.0
    mu = 0.0
    
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
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
    
    V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_k.data[:] = V
    
    print('--> pi_bubble')
    PI_wk = bubble_PI_wk(g0_wk)
    
    print('--> dynamical_screened_interaction_W')
    Wr_full_wk = dynamical_screened_interaction_W(PI_wk, V_k)

    Wr_dyn_wk, Wr_stat_k = split_into_dynamic_wk_and_constant_k(Wr_full_wk)

    Wr_dyn_wk_ref = Gf(mesh=Wr_full_wk.mesh, target_shape=[norb]*4)
    for w in Wr_full_wk.mesh.components[0]:
        iw = w.data_index
        Wr_dyn_wk_ref.data[iw,:] = Wr_full_wk.data[iw,:] - V_k.data[:]
  

    diff = Wr_dyn_wk.data[:] - Wr_dyn_wk_ref.data[:]
    print(np.max(np.abs(np.real(diff))))
    print(np.max(np.abs(np.imag(diff))))
    np.testing.assert_array_almost_equal(Wr_dyn_wk.data[:], Wr_dyn_wk_ref.data[:])
    
    diff = Wr_stat_k.data[:] - V_k.data[:]
    print(np.max(np.abs(np.real(diff))))
    print(np.max(np.abs(np.imag(diff))))
    np.testing.assert_array_almost_equal(Wr_stat_k.data[:], V_k.data[:])

def test_add_dynamic_and_static():
    print("== Add dynamic and static ==")
    
    nk = 3
    nb = 2
    wmin = -5.0
    wmax = 5.0
    beta = 2.0
    nw = 5
    lamb = 10.0
    eps = 1e-4

    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, nk])

    fmesh = MeshReFreq(wmin, wmax, nw)

    wmesh = MeshImFreq(beta, 'Fermion', nw)
    numesh = MeshImFreq(beta, 'Boson', nw)

    DLRwmesh = MeshDLRImFreq(beta, 'Fermion', lamb, eps)
    DLRnumesh = MeshDLRImFreq(beta, 'Boson', lamb, eps)

    def check_add_dynamic_static(freqmesh, kmesh, nb, rank, name=""):
        print("--> %s"%name)

        g_dyn_wk = Gf(mesh=MeshProduct(freqmesh, kmesh), target_shape=[nb]*rank)
        g_dyn_wk.data[:] = np.random.rand(*g_dyn_wk.data[:].shape)

        g_stat_k = Gf(mesh=kmesh, target_shape=[nb]*rank)
        g_stat_k.data[:] = np.random.rand(*g_stat_k.data[:].shape)

        g_wk = add_dynamic_and_static(g_dyn_wk, g_stat_k)

        g_ref_wk = Gf(mesh=MeshProduct(freqmesh, kmesh), target_shape=[nb]*rank)
        for w in freqmesh:
            wii = w.data_index
            g_ref_wk.data[wii,:] = g_dyn_wk.data[wii,:] + g_stat_k.data[:]

        np.testing.assert_array_almost_equal(g_wk.data[:], g_ref_wk.data[:])

    check_add_dynamic_static(fmesh, kmesh, nb, 2, "Real freq G")
    check_add_dynamic_static(fmesh, kmesh, nb, 4, "Real freq Chi")

    check_add_dynamic_static(wmesh, kmesh, nb, 2, "Linear Matsubara G")
    check_add_dynamic_static(numesh, kmesh, nb, 4, "Linear Matsubara Chi")

    check_add_dynamic_static(DLRwmesh, kmesh, nb, 2, "DLR Matsubara G")
    check_add_dynamic_static(DLRnumesh, kmesh, nb, 4, "DLR Matsubara Chi")

if __name__ == "__main__":
    test_split_into_dynamic_and_constant()
    test_add_dynamic_and_static()
    test_g_wk_to_g_mwk()
