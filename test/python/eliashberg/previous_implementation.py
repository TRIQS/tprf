# ----------------------------------------------------------------------

""" Goes through the steps of solving the linearized Eliashberg equation for singlet pairing in
RPA limit, saves the results and compares to previously established benchmark data. """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.lattice import eliashberg_product
from triqs_tprf.eliashberg import solve_eliashberg, allclose_by_scalar_multiplication

# ----------------------------------------------------------------------

from triqs_tprf.utilities import assert_parameter_collection_not_equal_model_parameters, write_TarGZ_HDFArchive, read_TarGZ_HDFArchive, show_version_info
import triqs_tprf.version as version

# ----------------------------------------------------------------------

def save_new_benchmarks(filename, p):
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    big_nw = 2*p.nw + 1
    eliashberg_ingredients_big = create_eliashberg_ingredients(p.alter(nw=big_nw))
    gamma_big = eliashberg_ingredients_big.gamma

    next_delta = eliashberg_product(gamma_big, g0_wk, g0_wk) 
    Es, eigen_modes = solve_eliashberg(gamma_big, g0_wk, product='SUM', solver='IRAM')

    p.next_delta = next_delta
    p.E = Es[0]
    p.eigen_mode = eigen_modes[0]

    write_TarGZ_HDFArchive(filename, p=p)

def test_next_delta(g0_wk, gamma_big, expected_next_delta):
    next_delta = eliashberg_product(gamma_big, g0_wk, g0_wk) 
    np.testing.assert_allclose(next_delta.data, expected_next_delta.data)

def test_solve_eliashberg(g0_wk, gamma_big, expected_E, expected_eigen_mode):
    Es, eigen_modes = solve_eliashberg(gamma_big, g0_wk, product='SUM', solver='IRAM')
    np.testing.assert_allclose(Es[0], expected_E) 
    assert allclose_by_scalar_multiplication(eigen_modes[0], expected_eigen_mode),\
                "Eigenvectors are not the same."


if __name__ == "__main__":
    p = ParameterCollection(
            benchmark_filename = "./eliashberg_benchmark.tar.gz",
            filename = 'eliashberg_benchmark_new.tar.gz',
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
    # For the eliashberg SUM procedure a Gamma with a twice as big w-mesh then the GF is needed.
    big_nw = 2*p.nw + 1
    eliashberg_ingredients_big = create_eliashberg_ingredients(p.alter(nw=big_nw))
    gamma_big = eliashberg_ingredients_big.gamma

    test_next_delta(g0_wk, gamma_big, p_benchmark.next_delta)
    test_solve_eliashberg(g0_wk, gamma_big, p_benchmark.E, p_benchmark.eigen_mode)

    print('\nThe benchmark data was obtained with %s.'%show_version_info(p_benchmark.version_info))
    print('\nThis (new) version with %s yields the same results!'%show_version_info(p.version_info))
