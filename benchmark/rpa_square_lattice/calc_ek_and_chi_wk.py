
# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq
from h5 import HDFArchive
from triqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from pyed.OperatorUtils import quartic_tensor_from_operator
from pyed.OperatorUtils import quartic_permutation_symmetrize

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    #n_k = (128, 128, 1) # need some parallelisation for this... FIXME
    #n_k = (96, 96, 1)
    n_k = (64, 64, 1)
    #n_k = (32, 32, 1)
    #n_k = (16, 16, 1)
    #n_k = (6, 6, 1)
    #n_k = (2, 2, 1)
    n_w = 40
    
    # ------------------------------------------------------------------
    # -- tight binding parameters

    beta = 10.0
    mu = 0.0
    t = 1.0
    U = 2.5 # ZERO !!
    
    h_loc = np.zeros((2, 2))        
    T = - t * np.eye(2)

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

    print('--> e_k')
    e_k = t_r.on_mesh_brillouin_zone(n_k)    

    print('--> g0_wk')
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_w)
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

    chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)

    # ------------------------------------------------------------------
    # -- RPA
    
    fundamental_operators = [c(0, 0), c(0, 1)]
    
    chi_wk_vec = []
    U_vec = np.arange(1.0, 5.0, 1.0)

    for u in U_vec:
        print('--> RPA chi_wk, U =', u)

        H_int = u * n(0, 0) * n(0, 1)
        U_int_abcd = quartic_tensor_from_operator(H_int, fundamental_operators)

        print(U_int_abcd.reshape((4, 4)).real)
        
        U_int_abcd = quartic_permutation_symmetrize(U_int_abcd)
        print(U_int_abcd.reshape((4, 4)).real)
        
        # -- Group in Gamma order cc^+cc^+ ( from c^+c^+cc )

        U_abcd = np.zeros_like(U_int_abcd)
        for a, b, c, d in itertools.product(list(range(U_abcd.shape[0])), repeat=4):
            U_abcd[a, b, c, d] = U_int_abcd[b, d, a, c]

        chi_wk = solve_rpa_PH(chi00_wk, U_abcd)

        chi_wk_vec.append(chi_wk)

    # -- Store to disk
    
    filename = 'data_ek_and_chi_wk.h5'
    
    with HDFArchive(filename, 'w') as arch:

        arch['n_k'] = n_k
        arch['e_k'] = e_k

        arch['U_vec'] = U_vec
        arch['chi_wk_vec'] = chi_wk_vec
        arch['chi00_wk'] = chi00_wk
