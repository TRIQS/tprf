from unittest.mock import patch, MagicMock

import numpy as np

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients

from triqs_tprf.eliashberg import solve_eliashberg


def test_no_initial_delta_input(g0_wk, gamma):
    with patch("triqs_tprf.eliashberg.semi_random_initial_delta") as patched:
        patched.return_value = g0_wk

        solve_eliashberg(gamma, g0_wk)

    patched.assert_called()


def test_initial_delta_input(g0_wk, gamma):
    with patch("triqs_tprf.eliashberg.semi_random_initial_delta") as patched:
        patched.return_value = g0_wk

        solve_eliashberg(gamma, g0_wk, initial_delta=g0_wk)

    patched.assert_not_called()


def test_tol_used_in_IRAM(g0_wk, gamma):
    expected_tol = 1e-4

    with patch("triqs_tprf.eliashberg.implicitly_restarted_arnoldi_method") as patched:
        patched.return_value = [0.0], [g0_wk.data.flatten()]

        solve_eliashberg(gamma, g0_wk, solver="IRAM", tol=expected_tol)

    kwargs = patched.call_args[-1]
    called_tol = kwargs["tol"]
    assert called_tol == expected_tol


def test_tol_used_in_PM(g0_wk, gamma):
    expected_tol = 1e-10

    with patch("triqs_tprf.eliashberg.power_method_LR") as patched:
        patched.return_value = 0.0, g0_wk.data.flatten()

        solve_eliashberg(gamma, g0_wk, solver="PM", tol=expected_tol)

    kwargs = patched.call_args[-1]
    called_tol = kwargs["tol"]
    assert called_tol == expected_tol


def test_call_eliashberg_product_fft(g0_wk, gamma):
    with patch("triqs_tprf.eliashberg.eliashberg_product_fft") as patched:
        patched.return_value = g0_wk

        solve_eliashberg(gamma, g0_wk, product="FFT")

    patched.assert_called()


def test_call_eliashberg_product_fft_constant(g0_wk, gamma):
    with patch("triqs_tprf.eliashberg.eliashberg_product_fft_constant") as patched:
        patched.return_value = g0_wk

        non_dynamic_gamma = 0 * gamma

        solve_eliashberg(non_dynamic_gamma, g0_wk, product="FFT")

    patched.assert_called()


def test_call_eliashberg_product(g0_wk, gamma):
    with patch("triqs_tprf.eliashberg.eliashberg_product") as patched:
        patched.return_value = g0_wk

        solve_eliashberg(gamma, g0_wk, product="SUM")

    patched.assert_called()


def test_wrong_input_for_product(g0_wk, gamma):
    wrong_input = "false"

    try:
        solve_eliashberg(gamma, g0_wk, product=wrong_input)
    except NotImplementedError as e:
        expected_message = (
            "There is no implementation of the eliashberg product called %s."
            % wrong_input
        )
        assert str(e) == expected_message


def test_call_IRAM_solver(g0_wk, gamma):
    with patch("triqs_tprf.eliashberg.implicitly_restarted_arnoldi_method") as patched:
        patched.return_value = [0.0], [g0_wk.data.flatten()]

        solve_eliashberg(gamma, g0_wk, solver="IRAM")

    patched.assert_called()


def test_call_PM_solver(g0_wk, gamma):
    with patch("triqs_tprf.eliashberg.power_method_LR") as patched:
        patched.return_value = 0.0, g0_wk.data.flatten()

        solve_eliashberg(gamma, g0_wk, solver="PM")

    patched.assert_called()


def test_wrong_input_for_solver(g0_wk, gamma):
    wrong_input = "false"

    try:
        solve_eliashberg(gamma, g0_wk, solver=wrong_input)
    except NotImplementedError as e:
        expected_message = "There is no solver called %s." % wrong_input
        assert str(e) == expected_message


def test_call_symmetrize_function(g0_wk, gamma):
    symmetrize_fct = MagicMock()
    symmetrize_fct.return_value = g0_wk

    solve_eliashberg(gamma, g0_wk, symmetrize_fct=symmetrize_fct)

    symmetrize_fct.assert_called()


def test_invalid_symmetrize_function(g0_wk, gamma):
    invalid_symmetrize_fct = lambda x: 1

    try:
        solve_eliashberg(gamma, g0_wk, symmetrize_fct=invalid_symmetrize_fct)
    except AttributeError as e:
        assert str(e) == "'int' object has no attribute 'data'"


def test_k_input(g0_wk, gamma):
    for k_input in [1, 3]:
        Es, evs = solve_eliashberg(gamma, g0_wk, k=k_input)
        assert len(Es) == k_input
        assert len(evs) == k_input


if __name__ == "__main__":
    p = ParameterCollection(
        dim=2,
        norb=1,
        t=1.0,
        mu=0.0,
        beta=1,
        U=1.0,
        Up=0.0,
        J=0.0,
        Jp=0.0,
        nk=2,
        nw=100,
    )

    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma

    test_no_initial_delta_input(g0_wk, gamma)
    test_initial_delta_input(g0_wk, gamma)
    test_tol_used_in_IRAM(g0_wk, gamma)
    test_tol_used_in_PM(g0_wk, gamma)
    test_call_eliashberg_product_fft(g0_wk, gamma)
    test_call_eliashberg_product_fft_constant(g0_wk, gamma)
    test_call_eliashberg_product(g0_wk, gamma)

    test_call_IRAM_solver(g0_wk, gamma)
    test_call_PM_solver(g0_wk, gamma)

    test_wrong_input_for_product(g0_wk, gamma)
    test_wrong_input_for_solver(g0_wk, gamma)
    test_call_symmetrize_function(g0_wk, gamma)
    test_invalid_symmetrize_function(g0_wk, gamma)
    test_k_input(g0_wk, gamma)
