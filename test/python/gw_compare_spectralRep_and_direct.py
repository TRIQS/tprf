# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk

from triqs_tprf.gw import gw_sigma
from triqs_tprf.gw import g0w_sigma

from triqs.gf import Gf, MeshImFreq
from triqs.gf.mesh_product import MeshProduct

# ----------------------------------------------------------------------

def test_gw_compare_direct_and_spectralRep():
    """ Compares the direct Matsubara sum method and the spectral-representation
    method of calculating the static GW self-energy 
    This is to test if the orbital index order is consistent between the two methods.
    The values used for the Coulomb interaction here are arbitrary and non-physical.
    Author: Yann in 't Veld (2023)"""

    nw = 1000
    norb = 2
    beta = 10.0
    A = 1.0 + 0.1j
    B = 0.9 - 0.05j
    C = 0.0
    C2 = 0.0
    D = 0.0
    mu = 0.0
    
    t = -1.0
    tp = -1.0
    t = np.array([[t,tp],[tp,t]])

    t_r = TBLattice(
        units = [(1, 0, 0)],
        hopping = {
            (+1,) : t,
            (-1,) : t,
            },
        orbital_positions = [(0,0,0)]*norb,
        )
    
    print('--> setup e_k')
    kmesh = t_r.get_kmesh(n_k=(1, 1, 1))
    e_k = t_r.fourier(kmesh)
    print(e_k.data[0,:])

    kmesh = e_k.mesh
    wmesh = MeshImFreq(beta, 'Fermion', nw)
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
    
    
    print('--> setup V_k')
    # Set V_k to some arbitrary values
    V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_k.data[:] = 0.0
    
    V_k.data[:,0,0,0,0] = A
    V_k.data[:,1,1,1,1] = A
    V_k.data[:,0,0,1,1] = B
    V_k.data[:,1,1,0,0] = np.conjugate(V_k.data[:,0,0,1,1])
    
    V_k.data[:,0,1,0,1] = C
    V_k.data[:,1,0,1,0] = np.conjugate(V_k.data[:,0,1,0,1])
    V_k.data[:,0,1,1,0] = C2
    V_k.data[:,1,0,0,1] = C2
    
    V_k.data[:,0,0,0,1] = D
    V_k.data[:,0,0,1,0] = D
    V_k.data[:,0,1,0,0] = D
    V_k.data[:,1,0,0,0] = D
    
    V_k.data[:,1,1,1,0] = D
    V_k.data[:,1,1,0,1] = D
    V_k.data[:,1,0,1,1] = D
    V_k.data[:,0,1,1,1] = D
    
    
    print('--> gw_sigma (direct)')
    sigma_k = gw_sigma(V_k, g0_wk)
    print(sigma_k.data[0,:])

    print('--> g0w_sigma (via spectral representation)')
    sigma_k_ref = g0w_sigma(mu=mu, beta=beta, e_k=e_k, v_k=V_k)
    print(sigma_k_ref.data[0,:])

    diff = sigma_k.data[:] - sigma_k_ref.data[:]
    print(np.max(np.abs(np.real(diff))))
    print(np.max(np.abs(np.imag(diff))))
    
    np.testing.assert_array_almost_equal(sigma_k.real.data[:], sigma_k_ref.real.data[:])
    np.testing.assert_array_almost_equal(sigma_k.imag.data[:], sigma_k_ref.imag.data[:])
    

if __name__ == "__main__":
    test_gw_compare_direct_and_spectralRep()

