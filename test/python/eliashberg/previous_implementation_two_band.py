
# ----------------------------------------------------------------------

""" Goes through the steps of solving the linearized Eliashberg equation for singlet pairing in
RPA limit in model with two orbitals, saves the results and compares to previously established 
benchmark data. """

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors
from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.lattice import eliashberg_product_fft
from triqs_tprf.eliashberg import preprocess_gamma_for_fft, solve_eliashberg
from triqs_tprf.eliashberg import allclose_by_scalar_multiplication

# ----------------------------------------------------------------------

from triqs_tprf.utilities import write_TarGZ_HDFArchive, read_TarGZ_HDFArchive, show_version_info
import triqs_tprf.version as version

# ----------------------------------------------------------------------

p = ParameterCollection(
        filename = 'eliashberg_benchmark_two_band_new.tar.gz',
        dim = 2,
        norbs = 2,
        t1 = 1.0,
        t2 = 0.5,
        t12 = 0.1,
        t21 = 0.1,
        mu = 0.0,
        beta = 1,
        U = 1.0,
        nk = 2,
        nw = 100,
        version_info = version.info,
        )

# -- Setup model, RPA susceptibilities and spin/charge interaction

full_units = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
all_nn_hoppings = list(itertools.product([-1, 0, 1], repeat=p.dim)) 
non_diagonal_hoppings = [ele for ele in all_nn_hoppings if sum(np.abs(ele)) == 1] 

# -- Create hopping matrix for two-band model
t = -np.array([[p.t1, p.t12], [p.t21, p.t2]])

H = TBLattice(
            units = full_units[:p.dim],
            hopping = {hop : t for hop in non_diagonal_hoppings},
            orbital_positions = [(0,0,0)]*p.norbs,
            )

e_k = H.on_mesh_brillouin_zone(n_k=[p.nk]*p.dim + [1]*(3-p.dim))

wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw)

g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)

chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw)

U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(p.norbs, p.U, 0, 0, 0)

chi_s = solve_rpa_PH(chi0_wk, U_s)
chi_c = solve_rpa_PH(chi0_wk, -U_c) # Minus for correct charge rpa equation

# -- The output of the following three functions shall be tested

gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)
Gamma_pp_dyn_tr, Gamma_pp_const_r = preprocess_gamma_for_fft(gamma) # This one is not tested
next_delta = eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g0_wk, g0_wk) 
E, eigen_modes = solve_eliashberg(gamma, g0_wk, product='FFT', solver='IRAM')

# -- Save results

p.gamma = gamma
p.next_delta = next_delta
p.E = E[0]
p.eigen_mode = eigen_modes[0]

write_TarGZ_HDFArchive(p.filename, p=p)

# -- Load benchmark data

filename = './eliashberg_benchmark_two_band.tar.gz'
p_benchmark = read_TarGZ_HDFArchive(filename)['p']

# -- Check if the benchmark data was calculated for the same model,
# -- otherwise a comparison does not make sense.

model_parameters = ['dim', 'norbs', 't1', 't2', 't12', 't21', 'mu', 'beta', 'U']

for model_parameter in model_parameters:
    try:
        run_time, benchmark = p[model_parameter], p_benchmark[model_parameter]
    except KeyError:
        raise AssertionError, "The model parameter %s does not exist."%model_parameter
    if (run_time != benchmark):
        error = 'The model of the benchmark and the one used now are not the same.\n' 
        error += '\t\tNow: {0} = {1}, benchmark: {0} = {2}.'.format(model_parameter, run_time,
                                                                                    benchmark)
        raise AssertionError, error

# -- Compare the results. Raise an error if the are not the same within a tolerance.

print('\nThe benchmark data was obtained with %s.'%show_version_info(p_benchmark.version_info))

np.testing.assert_allclose(p_benchmark.gamma.data, p.gamma.data)
np.testing.assert_allclose(p_benchmark.next_delta.data, p.next_delta.data)
np.testing.assert_allclose(p_benchmark.E, p.E) 
assert allclose_by_scalar_multiplication(p_benchmark.eigen_mode, p.eigen_mode),\
            "Eigenvectors are not the same."

print('\nThis (new) version with %s yields the same results!'%show_version_info(p.version_info))
