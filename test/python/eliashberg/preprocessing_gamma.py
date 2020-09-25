""" Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """
import numpy as np

from triqs.gf import MeshProduct, MeshImFreq, MeshBrillouinZone, MeshImTime, MeshCyclicLattice

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr
from triqs_tprf.eliashberg import preprocess_gamma_for_fft 

from triqs_tprf.utilities import assert_parameter_collection_not_equal_model_parameters, write_TarGZ_HDFArchive, read_TarGZ_HDFArchive
import triqs_tprf.version as version

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

def test_preprocess_gamma_for_fft_types(gamma):
    gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma) 

    assert type(gamma_dyn_tr.mesh) == MeshProduct
    assert type(gamma_dyn_tr.mesh[0]) == MeshImTime
    assert type(gamma_dyn_tr.mesh[1]) == MeshCyclicLattice

    assert type(gamma_const_r.mesh) == MeshCyclicLattice

def save_new_preprocess_gamma_for_fft_benchmark(filename, p):
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    gamma = eliashberg_ingredients.gamma

    gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma) 

    p.gamma = gamma
    p.gamma_dyn_tr = gamma_dyn_tr
    p.gamma_const_r = gamma_const_r 

    write_TarGZ_HDFArchive(filename, p=p)

def test_preprocess_gamma_for_fft_benchmark(gamma, p):
    p_benchmark = read_TarGZ_HDFArchive(p.benchmark_filename)['p']
    model_parameters_to_test = ['dim', 'norb', 't', 'mu', 'beta', 'U']
    assert_parameter_collection_not_equal_model_parameters(p, p_benchmark, model_parameters_to_test)

    np.testing.assert_allclose(gamma.data, p_benchmark.gamma.data, atol=1e-10)

    gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma) 
    np.testing.assert_allclose(gamma_dyn_tr.data, p_benchmark.gamma_dyn_tr.data, atol=1e-10)
    np.testing.assert_allclose(gamma_const_r.data, p_benchmark.gamma_const_r.data, atol=1e-10)


if __name__ == '__main__':
    p = ParameterCollection(
            filename = "preprocess_gamma_for_fft_benchmark_new.tar.gz",
            benchmark_filename = "preprocess_gamma_for_fft_benchmark.tar.gz",
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
            version_info = version.info,
            )
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    gamma = eliashberg_ingredients.gamma
    U_c = eliashberg_ingredients.U_c
    U_s = eliashberg_ingredients.U_s

    test_split_into_dynamic_wk_and_constant_k_mesh_types(gamma)
    test_split_into_dynamic_wk_and_constant_k_mesh_values(gamma, U_c, U_s)
    test_dynamic_and_constant_to_tr_mesh_types(gamma)
    test_preprocess_gamma_for_fft_types(gamma)

    #save_new_preprocess_gamma_for_fft_benchmark("preprocess_gamma_for_fft_benchmark.tar.gz", p)
    test_preprocess_gamma_for_fft_benchmark(gamma, p)
