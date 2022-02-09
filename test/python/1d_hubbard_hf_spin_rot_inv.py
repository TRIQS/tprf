# ----------------------------------------------------------------------

""" One dimensional Hubbard model solved with Hartree-Fock

Compare the result when applying a magnetic field in
many different directions to test rotational invariance.

Author: Hugo U.R. Strand (2018), hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.operators import n

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.hf_solver import HartreeSolver, HartreeFockSolver
from triqs_tprf.hf_response import HartreeResponse

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    beta = 1.3
    N_tot = 0.7
    n_k = (256, 1, 1)

    # -- One dimensional tight binding model

    t = 1.0
    h_loc = np.zeros((2, 2))
    T = - t * np.eye(2)
        
    t_r = TBLattice(
        units = [(1, 0, 0)],
        hopping = {
            # nearest neighbour hopping -t
            (0,): h_loc,
            (+1,): T,
            (-1,): T,
            },
        orbital_positions = [(0,0,0)] * 2,
        orbital_names = ['up', 'do'],
        )    

    kmesh = t_r.get_kmesh(n_k)
    e_k = t_r.fourier(kmesh)

    # -- Local double occupancy operator
    
    gf_struct = [[0, 2]]

    docc = n(0, 0) * n(0, 1)

    Sx = 0.5 * np.rot90(np.eye(2))
    Sy = 0.5 * np.rot90(np.diag([1.j, -1.j]))
    Sz = 0.5 * np.diag([1., -1.])

    print('Sx =\n', Sx)
    print('Sy =\n', Sy)
    print('Sz =\n', Sz)

    U = 2.4
    h = 0.5

    M0 = np.zeros((2, 2))
    mu0 = 0.5 * U # half-filling    
    H_int = U * docc

    e_k_z = e_k.copy()
    e_k_z.data[:] += h * Sz[None, ...]
    
    hs_z = HartreeSolver(e_k_z, beta, H_int=H_int, gf_struct=gf_struct)
    hs_z.solve_newton(N_target=N_tot, M0=M0, mu0=mu0)

    sz_exp = np.einsum('ab,ba', Sz, hs_z.density_matrix())
    print('sz_exp =', sz_exp)

    S_vec = [
        Sx, Sy, Sz,
        (Sx + Sy)/np.sqrt(2),
        (Sx + Sz)/np.sqrt(2),
        (Sy + Sz)/np.sqrt(2),
        (Sx + Sy)/np.sqrt(2),
        (Sx + Sy + Sz)/np.sqrt(3),
        ]
    
    for Si in S_vec:

        print('Si =\n', Si)
        
        e_k_h = e_k.copy()
        e_k_h.data[:] += h * Si[None, ...]
    
        hs = HartreeFockSolver(e_k_h, beta, H_int=H_int, gf_struct=gf_struct)
        hs.solve_newton(N_target=N_tot, M0=M0, mu0=mu0)
        
        si_exp = np.einsum('ab,ba', Si, hs.density_matrix())
        print('si_exp =', si_exp)

        np.testing.assert_almost_equal(si_exp, sz_exp)

        for attr in ['N_tot', 'E_int', 'E_kin', 'Omega0', 'E_tot', 'Omega']:
            v1 = getattr(hs_z, attr)
            v2 = getattr(hs, attr)
            print(attr, v1, v2) 
            np.testing.assert_almost_equal(v1, v2)
