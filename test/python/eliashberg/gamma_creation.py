""" Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """
import numpy as np

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs.gf import MeshProduct, MeshImFreq, MeshBrillouinZone

# from triqs_tprf.lattice import gamma_PP_spin_charge, gamma_PP_singlet, gamma_PP_triplet
from triqs_tprf.lattice import construct_phi_wk
from triqs_tprf.eliashberg import (
    construct_gamma_singlet_rpa,
    construct_gamma_triplet_rpa,
)


def test_phi_wk_mesh_type(chi_d, U_d):
    phi_d_wk = construct_phi_wk(chi_d, U_d)
    assert type(phi_d_wk.mesh) == MeshProduct
    assert type(phi_d_wk.mesh[0]) == MeshImFreq
    assert type(phi_d_wk.mesh[1]) == MeshBrillouinZone


def test_phi_wk_one_zero(chi_d, U_d):
    phi_d_wk = construct_phi_wk(chi_d, 0 * U_d)
    np.testing.assert_equal(phi_d_wk.data, 0.0)

    phi_d_wk = construct_phi_wk(0 * chi_d, U_d)
    np.testing.assert_equal(phi_d_wk.data, 0.0)


def test_gamma_singlet_mesh_type(chi_d, chi_m, U_d, U_m):
    phi_d_wk = construct_phi_wk(chi_d, U_d)
    phi_m_wk = construct_phi_wk(chi_m, U_m)

    gamma_singlet = construct_gamma_singlet_rpa(U_d, U_m, phi_d_wk, phi_m_wk)

    assert type(gamma_singlet.mesh) == MeshProduct
    assert type(gamma_singlet.mesh[0]) == MeshImFreq
    assert type(gamma_singlet.mesh[1]) == MeshBrillouinZone

def test_gamma_singlet_constant_only(chi_d, chi_m, U_d, U_m):
    phi_d_wk = construct_phi_wk(chi_d, U_d)
    phi_m_wk = construct_phi_wk(chi_m, U_m)

    gamma_singlet = construct_gamma_singlet_rpa(U_d, U_m, 0*phi_d_wk, 0*phi_m_wk)
    benchmark_value = 0.5 * U_d + 1.5 * U_m
    np.testing.assert_equal(gamma_singlet.data[0, 0], benchmark_value)


def test_gamma_triplet_mesh_type(chi_d, chi_m, U_d, U_m):
    phi_d_wk = construct_phi_wk(chi_d, U_d)
    phi_m_wk = construct_phi_wk(chi_m, U_m)

    gamma_triplet = construct_gamma_triplet_rpa(U_d, U_m, phi_d_wk, phi_m_wk)

    assert type(gamma_triplet.mesh) == MeshProduct
    assert type(gamma_triplet.mesh[0]) == MeshImFreq
    assert type(gamma_triplet.mesh[1]) == MeshBrillouinZone

def test_gamma_triplet_constant_only(chi_d, chi_m, U_d, U_m):
    phi_d_wk = construct_phi_wk(chi_d, U_d)
    phi_m_wk = construct_phi_wk(chi_m, U_m)

    gamma_triplet = construct_gamma_triplet_rpa(U_d, U_m, 0*phi_d_wk, 0*phi_m_wk)
    benchmark_value = -0.5 * U_d + 0.5 * U_m
    np.testing.assert_equal(gamma_triplet.data[0, 0], benchmark_value)

if __name__ == "__main__":
    p = ParameterCollection(
        dim=2,
        norb=2,
        t1=1.0,
        t2=0.5,
        t12=0.1,
        t21=0.1,
        mu=0.1,
        beta=1,
        U=1.0,
        Up=0.8,
        J=0.1,
        Jp=0.1,
        nk=3,
        nw=50,
    )

    eliashberg_ingredients = create_eliashberg_ingredients(p)
    chi_d = eliashberg_ingredients.chi_d
    chi_m = eliashberg_ingredients.chi_m
    U_d = eliashberg_ingredients.U_d
    U_m = eliashberg_ingredients.U_m

    test_phi_wk_mesh_type(chi_d, U_d)
    test_phi_wk_one_zero(chi_d, U_d)
    test_gamma_singlet_mesh_type(chi_d, chi_m, U_d, U_m)
    test_gamma_singlet_constant_only(chi_d, chi_m, U_d, U_m)
    test_gamma_triplet_mesh_type(chi_d, chi_m, U_d, U_m)
    test_gamma_triplet_constant_only(chi_d, chi_m, U_d, U_m)
