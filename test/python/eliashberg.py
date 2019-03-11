
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq, Gf, MeshProduct

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk

from triqs_tprf.lattice import eliashberg_product
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------
def test_eliashberg_product():
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    n_k = (2, 2, 1)
    nw = 10
    
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
    
    print '--> tight binding model'
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

    e_k = t_r.on_mesh_brillouin_zone(n_k)

    kmesh = e_k.mesh
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)

    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

    delta_in_wk = g0_wk.copy()
    delta_ref_wk = g0_wk.copy()

    gamma_pp = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=list(g0_wk.target_shape)*2)

    # -- Set gamma_pp and delta_ref_wk to something where we know the product

    gamma_pp.data[:] = np.random.random(gamma_pp.data.shape)
    delta_in_wk.data[:] = np.random.random(delta_in_wk.data.shape)
    delta_ref_wk.data[:] = np.random.random(delta_ref_wk.data.shape)

    # -- Compute the product
    
    print '--> Eliashberg product'
    delta_out_wk = eliashberg_product(gamma_pp, g0_wk, delta_in_wk)
    print 'done.'

    # -- Check that the result is the expected one

    #np.testing.assert_array_almost_equal(delta_out_wk.data, delta_ref_wk.data)

    E, eigen_modes = solve_eliashberg(gamma_pp, g0_wk)

    print 'E =', E

# ----------------------------------------------------------------------
if __name__ == '__main__':

    test_eliashberg_product()
