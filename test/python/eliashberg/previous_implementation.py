# ----------------------------------------------------------------------

""" Goes through the steps of solving the linearized Eliashberg equation for singlet pairing in
RPA limit, saves the results and compares to previously established benchmark data. """

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.lattice import eliashberg_product
from triqs_tprf.eliashberg import solve_eliashberg, allclose_by_scalar_multiplication

# ----------------------------------------------------------------------

from triqs_tprf.utilities import write_TarGZ_HDFArchive, read_TarGZ_HDFArchive, show_version_info
import triqs_tprf.version as version

# ----------------------------------------------------------------------

p = ParameterCollection(
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

eliashberg_ingredients = create_eliashberg_ingredients(p)
g0_wk = eliashberg_ingredients.g0_wk
gamma = eliashberg_ingredients.gamma
U_c = eliashberg_ingredients.U_c
U_s = eliashberg_ingredients.U_s

## A bigger w-mesh is needed to construct a Gamma with a twice as big w-mesh than GF
big_nw = 2*p.nw + 1
eliashberg_ingredients_big = create_eliashberg_ingredients(p.alter(nw=big_nw))
gamma_big = eliashberg_ingredients_big.gamma

# -- Setup model, RPA susceptibilities and spin/charge interaction

# -- The output of the following functions shall be tested
next_delta = eliashberg_product(gamma_big, g0_wk, g0_wk) 
E, eigen_modes = solve_eliashberg(gamma_big, g0_wk, product='SUM', solver='IRAM')

# -- Save results

p.gamma = gamma
p.next_delta = next_delta
p.E = E[0]
p.eigen_mode = eigen_modes[0]

write_TarGZ_HDFArchive(p.filename, p=p)

# -- Load benchmark data

filename = './eliashberg_benchmark.tar.gz'
p_benchmark = read_TarGZ_HDFArchive(filename)['p']

# -- Check if the benchmark data was calculated for the same model,
# -- otherwise a comparison does not make sense.

model_parameters = ['dim', 'norb', 't', 'mu', 'beta', 'U']

for model_parameter in model_parameters:
    run_time, benchmark = p[model_parameter], p_benchmark[model_parameter]
    if run_time != benchmark:
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
