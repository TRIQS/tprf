# ----------------------------------------------------------------------

""" Compare the implementations of the eliashberg products that use FFT.

The function `eliashberg_product_fft_constant` can only handle Gammas that are
constant in frequecny space, while `eliashberg_product_fft` can treat dynamic Gammas.
Here we test if both implementations give the same result for a Gamma that
is constant in momentum space.
This also tests the function 'split_into_dynamic_wk_and_constant_k', to
see if the split is done correctly.
"""

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import Gf, MeshImFreq, MeshProduct
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import eliashberg_product_fft, eliashberg_product_fft_constant
from triqs_tprf.eliashberg import semi_random_initial_delta, preprocess_gamma_for_fft

from triqs_tprf.tight_binding import create_model_for_tests
from triqs_tprf.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------

if __name__ == '__main__':

    p = ParameterCollection(
                            dim = 2,
                            norb = 1,
                            t = 2.0,
                            mu = 0.0,
                            beta = 5,
                            U = 1.0,
                            nk = 4,
                            nw = 200,
                                )

    H = create_model_for_tests(**p)
    e_k = H.on_mesh_brillouin_zone(n_k=(p.nk, p.nk, 1))

    wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw)
    g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)
        
    wmesh_boson = MeshImFreq(beta=p.beta, S='Boson', n_max=p.nw)
    gamma_pp_wk = Gf(mesh=MeshProduct(wmesh_boson, g0_wk.mesh[1]),
                  target_shape=g0_wk.target_shape*2)
    gamma_pp_wk.data[:] = np.random.rand(p.nk**2, 1, 1, 1, 1) 
    gamma_pp_dyn_tr, gamma_pp_const_r = preprocess_gamma_for_fft(gamma_pp_wk)

    initial_delta = semi_random_initial_delta(g0_wk)

    delta_1 = eliashberg_product_fft_constant(gamma_pp_const_r, g0_wk, initial_delta)
    delta_2 = eliashberg_product_fft(gamma_pp_dyn_tr, gamma_pp_const_r, g0_wk, initial_delta)

    np.testing.assert_allclose(delta_1.data, delta_2.data)

    print('The functions eliashberg_product_fft and eliashberg_product_fft_constant'
          ' yield the same result for a Gamma that is only constant in momentum space.'
          '\nThe function split_into_dynamic_wk_and_constant_k therefore also worked correcty.')
