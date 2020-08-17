import numpy as np

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients
from pytriqs.gf import MeshProduct, MeshImFreq, MeshBrillouinZone

from triqs_tprf.lattice import gamma_PP_spin_charge, gamma_PP_singlet, gamma_PP_triplet

def test_gamma_PP_spin_charge_mesh_types(chi_c, chi_s, U_c, U_s):
    gamma = gamma_PP_spin_charge(chi_c, chi_s, U_c, U_s, 0.0, 0.0)
    assert type(gamma.mesh) == MeshProduct
    assert type(gamma.mesh[0]) == MeshImFreq 
    assert type(gamma.mesh[1]) == MeshBrillouinZone 

def test_gamma_PP_spin_charge_zero_input(chi_c, chi_s, U_c, U_s):
    gamma = gamma_PP_spin_charge(0.0*chi_c, 0.0*chi_s, 0.0*U_c, 0.0*U_s, 0.0, 0.0)
    np.testing.assert_equal(gamma.data, 0)

def test_gamma_PP_spin_charge_only_constant(chi_c, chi_s, U_c, U_s):
    gamma = gamma_PP_spin_charge(chi_c, chi_s, U_c, U_s, 0.0, 0.0)
    np.testing.assert_equal(gamma.data[0, 0], 0.5*(U_c + U_s))

def test_gamma_PP_singlet_mesh_type(chi_c, chi_s, U_c, U_s):
    gamma_singlet = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)
    assert type(gamma_singlet.mesh) == MeshProduct
    assert type(gamma_singlet.mesh[0]) == MeshImFreq 
    assert type(gamma_singlet.mesh[1]) == MeshBrillouinZone 

def test_gamma_PP_singlet_value(chi_c, chi_s, U_c, U_s):
    gamma = gamma_PP_spin_charge(chi_c, chi_s, U_c, U_s, -1, 3)
    gamma_singlet = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)
    np.testing.assert_equal(gamma.data, gamma_singlet.data)

def test_gamma_PP_triplet_mesh_type(chi_c, chi_s, U_c, U_s):
    gamma_triplet = gamma_PP_triplet(chi_c, chi_s, U_c, U_s)
    assert type(gamma_triplet.mesh) == MeshProduct
    assert type(gamma_triplet.mesh[0]) == MeshImFreq 
    assert type(gamma_triplet.mesh[1]) == MeshBrillouinZone 

def test_gamma_PP_triplet_value(chi_c, chi_s, U_c, U_s):
    gamma = gamma_PP_spin_charge(chi_c, chi_s, U_c, U_s, -1, -1)
    gamma_triplet = gamma_PP_triplet(chi_c, chi_s, U_c, U_s)
    np.testing.assert_equal(gamma.data, gamma_triplet.data)

if __name__ == "__main__":
    p = ParameterCollection(
            dim = 2,
            norb = 2,
            t1 = 1.0,
            t2 = 0.5,
            t12 = 0.1,
            t21 = 0.1,
            mu = 0.1,
            beta = 1,
            U = 1.0,
            Up = 0.8,
            J = 0.1,
            Jp = 0.1,
            nk = 3,
            nw = 50,
            )

    eliashberg_ingredients = create_eliashberg_ingredients(p)
    chi_c = eliashberg_ingredients.chi_c
    chi_s = eliashberg_ingredients.chi_s
    U_c = eliashberg_ingredients.U_c
    U_s = eliashberg_ingredients.U_s

    test_gamma_PP_spin_charge_mesh_types(chi_c, chi_s, U_c, U_s)
    test_gamma_PP_spin_charge_zero_input(chi_c, chi_s, U_c, U_s)
    test_gamma_PP_spin_charge_only_constant(chi_c, chi_s, U_c, U_s)
    test_gamma_PP_singlet_mesh_type(chi_c, chi_s, U_c, U_s)
    test_gamma_PP_singlet_value(chi_c, chi_s, U_c, U_s)
    test_gamma_PP_triplet_mesh_type(chi_c, chi_s, U_c, U_s)
    test_gamma_PP_triplet_value(chi_c, chi_s, U_c, U_s)
