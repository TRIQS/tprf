
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq, MeshReFreq, Idx

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_fk

# ----------------------------------------------------------------------

from triqs_tprf.lattice import lindhard_chi00

# ----------------------------------------------------------------------
def test_square_lattice_chi00_realfreq():
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    n_k = (2, 2, 1)
    nw_g = 3
    wmin = -5.0
    wmax = 5.0
    nw = 10
    delta = 0.00001
    
    # ------------------------------------------------------------------
    # -- tight binding parameters

    beta = 40.0
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

    fmesh = MeshReFreq(wmin, wmax, nw_g)

    print('--> g0_fk')
    g0_fk = lattice_dyson_g0_fk(mu=mu, e_k=e_k, mesh=fmesh, delta=delta)

    # ------------------------------------------------------------------
    # -- chi00 in real frequencies at w=0
    
    print('--> chi00_fk')
    chi00_fk_analytic = lindhard_chi00(e_k=e_k, mesh=fmesh, beta=beta, mu=mu, delta=delta)

    f0_ind = len(fmesh)//2
    assert np.allclose( list(fmesh.values())[f0_ind], 0.0 )
    chi00_f0k = chi00_fk_analytic.data[f0_ind,:]
  
    # ------------------------------------------------------------------
    # -- chi00 in Matsubara frequencies at iw=0
  
    print('--> chi00_wk')
    wmesh_bose = MeshImFreq(beta, "Boson", nw)
    chi00_wk_analytic = lindhard_chi00(e_k=e_k, mesh=wmesh_bose, mu=mu)

    wmesh = chi00_wk_analytic.mesh.components[0]
    w0_ind = len(wmesh)//2
    assert np.allclose( list(wmesh.values())[w0_ind], 0.0 ) 
    chi00_w0k = chi00_wk_analytic.data[w0_ind,:]    
    
    # ------------------------------------------------------------------
    # -- compare at zero freq.
   
    np.testing.assert_array_almost_equal(
        chi00_f0k.data, chi00_w0k.data, decimal=5)

# ----------------------------------------------------------------------
if __name__ == '__main__':

    test_square_lattice_chi00_realfreq()
