import numpy as np

from pytriqs.gf import MeshProduct, MeshImFreq, MeshBrillouinZone, MeshImTime, MeshCyclicLattice

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.eliashberg import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr

def test_split_into_dynamic_wk_and_constant_k_mesh_types(gamma):
    gamma_dyn, gamma_const = split_into_dynamic_wk_and_constant_k(gamma)

    assert type(gamma_dyn.mesh) == MeshProduct
    assert type(gamma_dyn.mesh[0]) == MeshImFreq
    assert type(gamma_dyn.mesh[1]) == MeshBrillouinZone

    assert type(gamma_const.mesh) == MeshBrillouinZone

def test_split_into_dynamic_wk_and_constant_k_mesh_values(gamma, U_c, U_s):
    gamma_dyn, gamma_const = split_into_dynamic_wk_and_constant_k(gamma)

    analytical_constant_expression = 0.5*(U_s + U_c)
    np.testing.assert_allclose(gamma_const.data[0], analytical_constant_expression)

    gamma_without_constant_part = gamma.data - analytical_constant_expression 
    np.testing.assert_allclose(gamma_without_constant_part, gamma_dyn.data, atol=1e-12)

def test_dynamic_and_constant_to_tr_mesh_types(gamma):
    gamma_dyn, gamma_const = split_into_dynamic_wk_and_constant_k(gamma)
    gamma_dyn_tr, gamma_const_r = dynamic_and_constant_to_tr(gamma_dyn, gamma_const)

    assert type(gamma_dyn_tr.mesh) == MeshProduct
    assert type(gamma_dyn_tr.mesh[0]) == MeshImTime
    assert type(gamma_dyn_tr.mesh[1]) == MeshCyclicLattice

    assert type(gamma_const_r.mesh) == MeshCyclicLattice

if __name__ == '__main__':
    p = ParameterCollection(
            dim = 1,
            norb = 2,
            t = 2.0,
            mu = 0.0,
            beta = 5,
            U = 1.0,
            Up = 0.8,
            J = 0.1,
            Jp = 0.1,
            nk = 3,
            nw = 500,
            )
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    gamma = eliashberg_ingredients.gamma
    U_c = eliashberg_ingredients.U_c
    U_s = eliashberg_ingredients.U_s

    test_split_into_dynamic_wk_and_constant_k_mesh_types(gamma)
    test_split_into_dynamic_wk_and_constant_k_mesh_values(gamma, U_c, U_s)
    test_dynamic_and_constant_to_tr_mesh_types(gamma)
