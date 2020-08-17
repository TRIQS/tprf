from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients

from triqs_tprf.eliashberg import solve_eliashberg


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

def test_wrong_input_for_solver(g0_wk, gamma):
    wrong_input = "false"

    try:
        solve_eliashberg(gamma, g0_wk, solver=wrong_input)
    except NotImplementedError as e:
        expected_message = (
            "There is no solver called %s."
            % wrong_input
        )
        assert str(e) == expected_message

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

    test_wrong_input_for_product(g0_wk, gamma)
    test_wrong_input_for_solver(g0_wk, gamma)
