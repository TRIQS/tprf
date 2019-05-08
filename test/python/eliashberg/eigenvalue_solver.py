# ----------------------------------------------------------------------

""" Compare the summation implementation of the linearized Eliashberg product
and the one using Fourier transformations.
"""

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, MeshImFreq, Idx

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk, solve_rpa_PH
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------

def run_solve_eliashberg(p):

    # -- Setup model, RPA susceptibilities and spin/charge interaction

    full_units = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    all_nn_hoppings = list(itertools.product([-1, 0, 1], repeat=p.dim)) 
    non_diagonal_hoppings = [ele for ele in all_nn_hoppings if sum(np.abs(ele)) == 1] 

    t = -p.t * np.eye(p.norbs)

    H = TBLattice(
                units = full_units[:p.dim],
                hopping = {hop : t for hop in non_diagonal_hoppings},
                orbital_positions = [(0,0,0)]*p.norbs,
                )

    e_k = H.on_mesh_brillouin_zone(n_k=[p.nk]*p.dim + [1]*(3-p.dim))

    # A bigger w-mesh is needed to construct a Gamma with a twice as big w-mesh than GF

    wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw)
    wmesh_big = MeshImFreq(beta=p.beta, S='Fermion', n_max=int(p.big_factor*p.nw)+1)

    g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)
    g0_wk_big = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh_big)

    chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw)
    chi0_wk_big = imtime_bubble_chi0_wk(g0_wk_big, nw=int(p.big_factor*p.nw)+1)

    U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(p.norbs, p.U, p.Up, p.J,p.Jp)

    chi_s = solve_rpa_PH(chi0_wk, U_s)
    chi_c = solve_rpa_PH(chi0_wk, -U_c) # Minus for correct charge rpa equation
    chi_s_big = solve_rpa_PH(chi0_wk_big, U_s)
    chi_c_big = solve_rpa_PH(chi0_wk_big, -U_c) # Minus for correct charge rpa equation

    gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)
    gamma_big = gamma_PP_singlet(chi_c_big, chi_s_big, U_c, U_s)

    if p.product == 'SUM':
        gamma = gamma_big

    if p.fit_const:
        gamma_const = None
    else:
        gamma_const = 0.5*(U_s + U_c)
    
    Es, eigen_modes = solve_eliashberg(gamma, g0_wk, Gamma_pp_const_k=gamma_const, 
                                        product=p.product, solver=p.solver)
    
    return Es, eigen_modes

#================================================================================ 

if __name__ == '__main__':

    p = ParameterCollection(
            dim = 1,
            norbs = 1,
            t = 1.0,
            mu = 0.0,
            beta = 5,
            U = 1.0,
            Up = 0.8,
            J = 0.1,
            Jp = 0.1,
            nk = 4,
            nw = 350,
            fit_const = False,
            big_factor = 2,
            product = 'FFT',
            solver = 'PM',
            )

    Es_pm, eigen_modes_pm = run_solve_eliashberg(p)

    Es_iram, eigen_modes_iram = run_solve_eliashberg(p.copy(solver='IRAM'))

    print(Es_pm[0], Es_iram[0])

    np.testing.assert_allclose(Es_pm[0], Es_iram[0])

    try:
        np.testing.assert_allclose(eigen_modes_pm[0].data, eigen_modes_iram[0].data, atol=1e-8)
    except AssertionError:
        np.testing.assert_allclose(-eigen_modes_pm[0].data, eigen_modes_iram[0].data, atol=1e-8)

    print('Both solvers yield the same results.')

