
# ----------------------------------------------------------------------

""" Goes through the steps of solving the linearized Eliashberg equation for singlet pairing in
RPA limit in model with two orbitals, saves the results and compares to previously established 
benchmark data.

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq

from triqs_tprf.ParameterCollection import ParameterCollection

from triqs_tprf.tight_binding import create_model_for_tests

from triqs_tprf.lattice import eliashberg_product_fft
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.eliashberg import preprocess_gamma_for_fft, solve_eliashberg
from triqs_tprf.eliashberg import allclose_by_scalar_multiplication

# ----------------------------------------------------------------------

from triqs_tprf.utilities import assert_parameter_collection_not_equal_model_parameters, write_TarGZ_HDFArchive, read_TarGZ_HDFArchive, show_version_info
import triqs_tprf.version as version

# ----------------------------------------------------------------------

def save_new_benchmarks(filename, p):
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma

    Gamma_pp_dyn_tr, Gamma_pp_const_r = preprocess_gamma_for_fft(gamma)
    next_delta = eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g0_wk, g0_wk) 
    Es, eigen_modes = solve_eliashberg(gamma, g0_wk, product='FFT', solver='IRAM')

    p.next_delta = next_delta
    p.E = Es[0]
    p.eigen_mode = eigen_modes[0]

    write_TarGZ_HDFArchive(filename, p=p)

def test_next_delta(g0_wk, gamma, expected_next_delta):
    Gamma_pp_dyn_tr, Gamma_pp_const_r = preprocess_gamma_for_fft(gamma)
    next_delta = eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g0_wk, g0_wk) 
    np.testing.assert_allclose(next_delta.data, expected_next_delta.data, atol=10e-12)

def test_solve_eliashberg(g0_wk, gamma, expected_E, expected_eigen_mode):
    Es, eigen_modes = solve_eliashberg(gamma, g0_wk, product='FFT', solver='IRAM')
    np.testing.assert_allclose(Es[0], expected_E) 
    assert allclose_by_scalar_multiplication(eigen_modes[0], expected_eigen_mode),\
                "Eigenvectors are not the same."

if __name__ == "__main__":
    p = ParameterCollection(
            benchmark_filename = "./eliashberg_benchmark_one_band.tar.gz",
            filename = 'eliashberg_benchmark_one_band_new.tar.gz',
            dim = 2,
            norb = 1,
            t = 1.0,
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

    #save_new_benchmarks(p.filename, p)

    p_benchmark = read_TarGZ_HDFArchive(p.benchmark_filename)['p']
    model_parameters_to_test = ['dim', 'norb', 't', 'mu', 'beta', 'U']
    assert_parameter_collection_not_equal_model_parameters(p, p_benchmark, model_parameters_to_test)

    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma

    test_next_delta(g0_wk, gamma, p_benchmark.next_delta)
    test_solve_eliashberg(g0_wk, gamma, p_benchmark.E, p_benchmark.eigen_mode)

    print(('\nThe benchmark data was obtained with %s.'%show_version_info(p_benchmark.version_info)))
    print(('\nThis (new) version with %s yields the same results!'%show_version_info(p.version_info)))
