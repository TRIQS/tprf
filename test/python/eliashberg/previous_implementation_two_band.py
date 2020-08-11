
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

from triqs_tprf.tight_binding import create_model_for_tests

from triqs_tprf.lattice import eliashberg_product_fft
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.eliashberg import preprocess_gamma_for_fft, solve_eliashberg
from triqs_tprf.eliashberg import allclose_by_scalar_multiplication

# ----------------------------------------------------------------------

from triqs_tprf.utilities import write_TarGZ_HDFArchive, read_TarGZ_HDFArchive, show_version_info
import triqs_tprf.version as version

# ----------------------------------------------------------------------

p = ParameterCollection(
        filename = 'eliashberg_benchmark_two_band_new.tar.gz',
        dim = 2,
        norb = 2,
        t1 = 1.0,
        t2 = 0.5,
        t12 = 0.1,
        t21 = 0.1,
        mu = 0.0,
        beta = 1,
        U = 1.0,
        Up = 0.0,
        J = 0.0,
        Jp = 0.0,
        nk = 2,
        nw = 100,
        version_info = version.info,
        )

eliashberg_ingredients = create_eliashberg_ingredients(p)
g0_wk = eliashberg_ingredients.g0_wk
gamma = eliashberg_ingredients.gamma

Gamma_pp_dyn_tr, Gamma_pp_const_r = preprocess_gamma_for_fft(gamma)
# -- The output of the following functions shall be tested
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

model_parameters = ['dim', 'norb', 't1', 't2', 't12', 't21', 'mu', 'beta', 'U']

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

np.testing.assert_allclose(p_benchmark.gamma.data, p.gamma.data, atol=1e-9)
np.testing.assert_allclose(p_benchmark.next_delta.data, p.next_delta.data, atol=1e-9)
np.testing.assert_allclose(p_benchmark.E, p.E) 
assert allclose_by_scalar_multiplication(p_benchmark.eigen_mode, p.eigen_mode),\
            "Eigenvectors are not the same."

print('\nThis (new) version with %s yields the same results!'%show_version_info(p.version_info))
