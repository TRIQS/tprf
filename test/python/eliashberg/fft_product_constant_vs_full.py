# ----------------------------------------------------------------------

""" Compare the implementations of the eliashberg products that use FFT.

The function `eliashberg_product_fft_constant` can only handle Gammas that are
constant in frequecny space, while `eliashberg_product_fft` can treat dynamic Gammas.
Here we test if both implementations give the same result for a Gamma that
is constant in momentum space.
This also tests the function 'split_into_dynamic_wk_and_constant_k', to
see if the split is done correctly.

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import Gf, MeshImFreq, MeshProduct
from triqs_tprf.lattice import eliashberg_product_fft, eliashberg_product_fft_constant
from triqs_tprf.eliashberg import semi_random_initial_delta, preprocess_gamma_for_fft
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------

def test_eliashberg_product_fft_constant(g0_wk, gamma):
    gamma.data[:] = np.random.rand(*gamma.data.shape[1:])
    gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma)

    initial_delta = semi_random_initial_delta(g0_wk)

    delta_1 = eliashberg_product_fft_constant(gamma_const_r, g0_wk, initial_delta)
    delta_2 = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, initial_delta)

    np.testing.assert_allclose(delta_1.data, delta_2.data)

    print('The functions eliashberg_product_fft and eliashberg_product_fft_constant'
          ' yield the same result for a Gamma that is only constant in momentum space.'
          '\nThe function split_into_dynamic_wk_and_constant_k therefore also worked correcty.')

if __name__ == '__main__':

    p = ParameterCollection(
                            dim = 2,
                            norb = 1,
                            t = 2.0,
                            mu = 0.0,
                            beta = 5,
                            U = 1.0,
                            Up = 0.0,
                            J = 0.0,
                            Jp = 0.0,
                            nk = 4,
                            nw = 200,
                                )
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma

    test_eliashberg_product_fft_constant(g0_wk, gamma)
