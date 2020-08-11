# ----------------------------------------------------------------------

""" Compare the output of the implemented eigenvalue solver:
The Power Method and the Implicitly Restarted Arnoldi Method.

"""

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, MeshImFreq, Idx
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.eliashberg import solve_eliashberg, semi_random_initial_delta
from triqs_tprf.eliashberg import allclose_by_scalar_multiplication

# ----------------------------------------------------------------------

def run_solve_eliashberg(p):
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma

    initial_delta = semi_random_initial_delta(g0_wk, seed=1337)
    Es, eigen_modes = solve_eliashberg(gamma, g0_wk, 
                                        product="FFT", solver=p.solver, initial_delta=initial_delta)
    
    return Es, eigen_modes

#================================================================================ 

if __name__ == '__main__':

    p = ParameterCollection(
            dim = 1,
            norb = 1,
            t = 1.0,
            mu = 0.0,
            beta = 5,
            U = 1.0,
            Up = 0.8,
            J = 0.1,
            Jp = 0.1,
            nk = 4,
            nw = 200,
            solver = 'PM',
            )

    Es_pm, eigen_modes_pm = run_solve_eliashberg(p)

    Es_iram, eigen_modes_iram = run_solve_eliashberg(p.alter(solver='IRAM'))

    print(Es_pm[0], Es_iram[0])
    np.testing.assert_allclose(Es_pm[0], Es_iram[0])

    assert allclose_by_scalar_multiplication(eigen_modes_pm[0], eigen_modes_iram[0]),\
            "Eigenvectors are not the same."

    print('Both solvers yield the same results.')
