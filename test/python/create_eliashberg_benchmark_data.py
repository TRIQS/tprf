# ----------------------------------------------------------------------

""" Create benchmark file for the linearized Eliashberg equation

Goes through the steps of solving the linearized Eliashberg equation for singlet pairing in
RPA limit and saves the input and output to compare with succeding implementations. """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq
from pytriqs.archive import HDFArchive

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors
from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.lattice import eliashberg_product
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------

from utilities import get_git_revision_short_hash, write_TarGZ_HDFArchive

def save_eliashberg_benchmark_data(p):
    """ 
    Parameters:

    p : ParameterCollection, see below to see what it should contain
    """

    # -- Setup model, RPA susceptibilities and spin/charge interaction

    H = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        hopping = {
            # nearest neighbour hopping -t
            ( 0,+1): -p.t * np.eye(p.norbs),
            ( 0,-1): -p.t * np.eye(p.norbs),
            (+1, 0): -p.t * np.eye(p.norbs),
            (-1, 0): -p.t * np.eye(p.norbs),
            },
        orbital_positions = [(0,0,0)]*p.norbs,
        )

    e_k = H.on_mesh_brillouin_zone(n_k = (p.nk, p.nk, 1))

    wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw_gf)
    g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)

    chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw_chi0)

    U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(p.norbs, p.U, 0, 0, 0)

    chi_s = solve_rpa_PH(chi00_wk, U_s)
    chi_c = solve_rpa_PH(chi00_wk, -U_c) # Minus for correct charge rpa equation

    # -- The output of the following three functions shall be tested

    gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)
    next_delta = eliashberg_product(gamma, g0_wk, g0_wk) 
    E, eigen_modes = solve_eliashberg(gamma, g0_wk)
    
    # -- Input
    p.e_k = e_k
    p.g0_wk = g0_wk
    p.chi_c = chi_c
    p.chi_s = chi_s
    p.U_c = U_c
    p.U_s = U_s
    # -- Output
    p.gamma = gamma
    p.next_delta = next_delta
    p.E = E[0]
    p.eigen_mode = eigen_modes[0]

    write_TarGZ_HDFArchive(p.filename, p=p)

# ----------------------------------------------------------------------

if __name__ == '__main__':

    p = ParameterCollection(
            filename = 'eliashberg_benchmark.tar.gz',
            norbs = 1,
            t = 1.0,
            mu = 0.0,
            beta = 1,
            U = 1.0,
            nk = 2,
            nw_gf = 100,
            nw_chi0 = 100,
            tprf_hash = get_git_revision_short_hash()
            )

    save_eliashberg_benchmark_data(p)
