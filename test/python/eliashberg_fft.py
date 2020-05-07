# ----------------------------------------------------------------------

""" Compare the naive implementation of the linearized Eliashberg product
and the one using Fourier transformations.
This test is quite computational intensive, because of the inefficiency of
the naive implementations. In the future the Fourier transformation 
implementation will subsitute the naive one.
"""

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs.gf import Gf, MeshImFreq

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk, solve_rpa_PH
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.lattice import eliashberg_product, eliashberg_product_fft
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr
from triqs_tprf.eliashberg import solve_eliashberg, solve_eliashberg_fft
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors

# ----------------------------------------------------------------------

p = ParameterCollection(
        dim = 1,
        norbs = 1,
        t = 1.0,
        mu = 0.0,
        beta = 5,
        U = 1.0,
        nk = 4,
        nw = 500,
        )

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
big_factor = 2.0

wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw)
wmesh_big = MeshImFreq(beta=p.beta, S='Fermion', n_max=int(big_factor*p.nw))

g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)
g0_wk_big = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh_big)

chi0_wk_big = imtime_bubble_chi0_wk(g0_wk_big, nw=int(big_factor*p.nw)+1)

U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(p.norbs, p.U, 0, 0, 0)

chi_s_big = solve_rpa_PH(chi0_wk_big, U_s)
chi_c_big = solve_rpa_PH(chi0_wk_big, -U_c) # Minus for correct charge rpa equation

gamma_big = gamma_PP_singlet(chi_c_big, chi_s_big, U_c, U_s)

# -- Preprocess gamma for the FFT implementation

gamma_dyn_wk, gamma_const_k = split_into_dynamic_wk_and_constant_k(gamma_big)
gamma_dyn_tr, gamma_const_r = dynamic_and_constant_to_tr(gamma_dyn_wk, gamma_const_k)

# -- Test the Eliashberg equation

next_delta = eliashberg_product(gamma_big, g0_wk, g0_wk)
next_delta_fft = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, g0_wk)

np.testing.assert_allclose(next_delta.data, next_delta_fft.data, atol=1e-7)

Es, eigen_modes = solve_eliashberg(gamma_big, g0_wk)
Es_fft, eigen_modes_fft = solve_eliashberg_fft(gamma_big, g0_wk)

E = Es[0]
eigen_mode = eigen_modes[0]
E_fft = Es_fft[0]
eigen_mode_fft = eigen_modes_fft[0]

np.testing.assert_allclose(E, E_fft, atol=1e-7) 

try:
    np.testing.assert_allclose(eigen_mode.data, eigen_mode_fft.data, atol=1e-7) 
except AssertionError:
    np.testing.assert_allclose(-eigen_mode.data, eigen_mode_fft.data, atol=1e-7) 

print('\nSame results for both implementations of the linearized Eliashberg equation.')
