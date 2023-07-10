
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq, Idx, Gf
from triqs.gf.gf_factories import make_gf_dlr

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.gf.meshes import MeshDLRImFreq
from triqs_tprf.lattice import dlr_on_imfreq
from triqs_tprf.lattice import lindhard_chi00

# ----------------------------------------------------------------------


def compare_g_Dwk_and_g_wk(g_Dwk, g_wk, decimal=7):
    wmesh, kmesh = g_wk.mesh.components
    DLRwmesh = g_Dwk.mesh.components[0]

    g_ref_wk =Gf(mesh=g_wk.mesh, target_shape=g_wk.target_shape)
    for k in kmesh:
        g_Dw = Gf(mesh=DLRwmesh, target_shape=g_Dwk.target_shape)
        g_Dw.data[:] = g_Dwk.data[:,k.data_index,:]
        g_Dc = make_gf_dlr(g_Dw)
        g_ref_wk[:,k] = dlr_on_imfreq(g_Dc, wmesh)
    
    np.testing.assert_array_almost_equal(g_wk.data[:], g_ref_wk.data[:], decimal=decimal)

def test_square_lattice_chi00_dlr():
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    n_k = (2, 2, 1)
    nw = 50
    lamb = 10.
    eps = 1e-8

    # ------------------------------------------------------------------
    # -- tight binding parameters

    beta = 20.0
    mu = 0.0
    t = 1.0
    
    h_loc = np.array([
        [-0.3, -0.5],
        [-0.5, .4],
        ])
        
    T = - t * np.array([
        [1., 0.23],
        [0.23, 0.5],
        ])

    # ------------------------------------------------------------------
    # -- tight binding
    
    print('--> tight binding model')
    t_r = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        hopping = {
            # nearest neighbour hopping -t
            ( 0, 0): h_loc,
            ( 0,+1): T,
            ( 0,-1): T,
            (+1, 0): T,
            (-1, 0): T,
            },
        orbital_positions = [(0,0,0)]*2,
        orbital_names = ['up_0', 'do_0'],
        )

    kmesh = t_r.get_kmesh(n_k)
    e_k = t_r.fourier(kmesh)

    wmesh_bose = MeshImFreq(beta=beta, S='Boson', n_max=nw)
    DLRwmesh_bose = MeshDLRImFreq(beta, 'Boson', lamb, eps)

    
    print('--> chi00_wk analytic')
    chi00_wk_analytic = lindhard_chi00(e_k=e_k, mesh=wmesh_bose, mu=mu)
    chi00_Dwk_analytic = lindhard_chi00(e_k=e_k, mesh=DLRwmesh_bose, mu=mu)

    print('--> compare')
    compare_g_Dwk_and_g_wk(chi00_Dwk_analytic, chi00_wk_analytic)

# ----------------------------------------------------------------------
if __name__ == '__main__':
    test_square_lattice_chi00_dlr()

