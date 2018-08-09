
# ----------------------------------------------------------------------

import time
import itertools
import numpy as np

from scipy.linalg import expm

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive
from pytriqs.operators import n, c, c_dag, Operator, dagger
from pytriqs.statistics.histograms import Histogram

from pytriqs.operators.util.op_struct import set_operator_structure
from pytriqs.operators.util.U_matrix import U_matrix_kanamori, U_matrix
from pytriqs.operators.util.hamiltonians import h_int_kanamori

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq, Idx

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.lattice import solve_rpa_PH

# ----------------------------------------------------------------------

from triqs_tprf.hartree_solver import *

# ----------------------------------------------------------------------
if __name__ == '__main__':

    mu = 0.0    
    beta, n_w, n_k = 1.0, 100, 16

    U, J = 2.3, 0.4
    print 'U, J =', U, J

    spin_names = ['up', 'do']
    orb_names = [0, 1, 2]
    #orb_names = [0]
    norb = 2*len(orb_names)

    gf_struct = set_operator_structure(spin_names, orb_names, False) # orbital diag    
    U_mat, UPrime_mat = U_matrix_kanamori(n_orb=len(orb_names), U_int=U, J_hund=J)
    H_int = h_int_kanamori(
        spin_names, orb_names, U_mat, UPrime_mat, J_hund=J,
        off_diag=False, map_operator_structure=None, H_dump=None) # orbital diag

    # ------------------------------------------------------------------
    # -- Tightbinding model

    t = 1.0
    h_loc = np.kron(np.eye(2), np.diag([0., -0.1, 0.1]))
    T = - t * np.kron(np.eye(2), np.diag([0.01, 1., 1.]))

    t_r = TBLattice(
        units = [(1, 0, 0)],
        hopping = {
            # nearest neighbour hopping -t
            (0,): h_loc,
            (+1,): T,
            (-1,): T,
            },
        orbital_positions = [(0,0,0)] * norb,
        orbital_names = ['up_0', 'up_1', 'up_2', 'do_0', 'do_1', 'do_2'],
        )

    Sz = np.kron(np.diag([+0.5, -0.5]), np.eye(norb/2))

    if False:
        h_loc = np.kron(np.eye(2), np.diag([0.]))
        T = - t * np.kron(np.eye(2), np.diag([1.0]))

        t_r = TBLattice(
            units = [(1, 0, 0)],
            hopping = {
                # nearest neighbour hopping -t
                (0,): h_loc,
                (+1,): T,
                (-1,): T,
                },
            orbital_positions = [(0,0,0)] * norb,
            orbital_names = ['up_0', 'do_0'],
            )

        Sz = np.kron(np.diag([+0.5, -0.5]), np.eye((1)))

    print 'h_loc =\n', h_loc
    print 'T = \n', T
    print 'Sz =\n', Sz
    
    n_k = tuple([n_k, 1, 1])
    e_k = t_r.on_mesh_brillouin_zone(n_k)
    
    # ------------------------------------------------------------------
    # -- Hartree solver

    hs = DensityDensityMeanFieldSolver(e_k, H_int, gf_struct)
    hs.beta = beta
    hs.mu = 0.

    eps = 1e-9

    hr = HartreeResponse(hs, eps=eps)

    hr.chi0_SzSz = hr.bare_response(Sz, Sz)
    hr.chi_SzSz = hr.response(Sz, Sz)

    hfr = HartreeFockResponse(hs, eps=eps)

    hfr.chi0_SzSz = hfr.bare_response(Sz, Sz)
    hfr.chi_SzSz = hfr.response(Sz, Sz)

    np.testing.assert_almost_equal(hr.chi0_SzSz, hfr.chi0_SzSz)
    np.testing.assert_almost_equal(hr.chi_SzSz, hfr.chi_SzSz)
    
    # ------------------------------------------------------------------
    # -- Call TPRF chi0_wk bubble calc

    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_w)
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

    chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)

    chi0_abcd = chi0_wk[Idx(0), Idx(0,0,0)].copy()
    chi0_ab = extract_dens_dens(chi0_abcd)
    chi0_SzSz = np.einsum('ab,abcd,cd->', Sz, chi0_abcd, Sz)

    print 'chi0_abcd diff =', np.max(np.abs(chi0_abcd - hfr.chi0_abcd))
    np.testing.assert_almost_equal(chi0_abcd, hfr.chi0_abcd, decimal=6)
    
    print 'chi0_ab diff =', np.max(np.abs(chi0_ab - hr.chi0_ab))
    np.testing.assert_almost_equal(chi0_ab, hr.chi0_ab)
    
    np.testing.assert_almost_equal(chi0_SzSz, hr.chi0_SzSz)
        
    # ------------------------------------------------------------------

    chi_wk = solve_rpa_PH(chi0_wk, hs.U_abcd)

    chi_abcd = chi_wk[Idx(0), Idx(0,0,0)].copy()
    chi_ab = extract_dens_dens(chi_abcd)

    chi_SzSz = np.einsum('aa,ab,bb->', Sz, chi_ab, Sz)
    chi_SzSz_2 = np.einsum('ab,abcd,cd->', Sz, chi_abcd, Sz)

    print 'chi_abcd diff =', np.max(np.abs(chi_abcd - hfr.chi_abcd))
    np.testing.assert_almost_equal(chi_abcd, hfr.chi_abcd, decimal=6)

    print 'chi_ab diff =', np.max(np.abs(chi_ab - hr.chi_ab))
    np.testing.assert_almost_equal(chi_ab, hr.chi_ab)
    
    np.testing.assert_almost_equal(chi_SzSz, hr.chi_SzSz)
    np.testing.assert_almost_equal(chi_SzSz, hr.chi_SzSz)
    
