# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k

from triqs_tprf.gw import bubble_PI_wk
from triqs_tprf.gw import dynamical_screened_interaction_W

from triqs.gf import Gf, MeshImFreq, Idx, MeshImTime
from triqs.gf.mesh_product import MeshProduct

# ----------------------------------------------------------------------

def test_split_into_dynamic_and_constant():
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
        iw = w.linear_index
        Wr_dyn_wk_ref.data[iw,:] = Wr_full_wk.data[iw,:] - V_k.data[:]
  

    diff = Wr_dyn_wk.data[:] - Wr_dyn_wk_ref.data[:]
    print(np.max(np.abs(np.real(diff))))
    print(np.max(np.abs(np.imag(diff))))
    np.testing.assert_array_almost_equal(Wr_dyn_wk.data[:], Wr_dyn_wk_ref.data[:])
    
    diff = Wr_stat_k.data[:] - V_k.data[:]
    print(np.max(np.abs(np.real(diff))))
    print(np.max(np.abs(np.imag(diff))))
    np.testing.assert_array_almost_equal(Wr_stat_k.data[:], V_k.data[:])


if __name__ == "__main__":
    test_split_into_dynamic_and_constant()

