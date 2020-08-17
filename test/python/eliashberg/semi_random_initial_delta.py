import numpy as np

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients

from triqs_tprf.eliashberg import semi_random_initial_delta

def test_same_output_for_same_seed(g0_wk):
    random_delta1 = semi_random_initial_delta(g0_wk, seed=1)
    random_delta2 = semi_random_initial_delta(g0_wk, seed=1)

    np.testing.assert_equal(random_delta1.data,  random_delta2.data)

def test_different_output_for_different_seed(g0_wk):
    random_delta1 = semi_random_initial_delta(g0_wk, seed=1)
    random_delta2 = semi_random_initial_delta(g0_wk, seed=42)

    try:
        np.testing.assert_equal(random_delta1.data,  random_delta2.data)
        raise ValueError
    except AssertionError:
        pass

if __name__ == "__main__":

    p = ParameterCollection(
            dim = 2,
            norb = 2,
            t1 = 1.0,
            t2 = 0.5,
            t12 = 0.1,
            t21 = 0.1,
            mu = 0.0,
            beta = 1,
            U = 0.0,
            Up = 0.0,
            J = 0.0,
            Jp = 0.0,
            nk = 2,
            nw = 100,
            )

    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk

    test_same_output_for_same_seed(g0_wk)
    test_different_output_for_different_seed(g0_wk)
