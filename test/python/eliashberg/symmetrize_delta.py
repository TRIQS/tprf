# ----------------------------------------------------------------------

""" Test the symmetrizing feature of the solve_eliashberg function

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """

# ----------------------------------------------------------------------

import itertools
import functools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import Idx

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------

from triqs_tprf.symmetries import enforce_symmetry, check_symmetry

# ----------------------------------------------------------------------

def test_symmetry_constraint_of_solve_eliashberg(g0_wk, gamma, do_plot=False):
    variables = ["frequency", "momentum", "orbital"]
    all_symmetries = [('even', 'odd', 'even'), ('odd', 'even', 'even')]

    for symmetries in all_symmetries:
        symmetrize_fct = functools.partial(enforce_symmetry, 
                                       variables=variables,
                                       symmetries=symmetries)

        E, eigen_modes = solve_eliashberg(gamma, g0_wk, product='FFT', solver='IRAM',
                                            symmetrize_fct=symmetrize_fct)

        translate_symmetries = {"even" : +1, "odd" : -1, None : None}
        expected_symmetries = {variable : translate_symmetries[symmetry] \
                        for (variable, symmetry) in zip(variables, symmetries)}

        for delta in eigen_modes:
            produced_symmetries = check_symmetry(delta)
            if not expected_symmetries == produced_symmetries:
                raise AssertionError("Incorrect symmetries were produced.")

            if do_plot:
                plot_delta(delta)

def plot_delta(delta):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 3)

    vmax = np.max(np.abs(delta[Idx(0),:].data))

    for orb1, orb2 in itertools.product(list(range(p.norb)), repeat=2):
        shape = (p.nk, p.nk, p.norb, p.norb)
        data = delta[Idx(0), :].data.reshape(shape)
        plt.sca(axes[orb1,orb2])
        plt.imshow(data[:,:,orb1,orb2].real, cmap="RdBu_r",
                                                 vmax=vmax, vmin=-vmax)
        plt.colorbar()

    plt.sca(axes[-1,-1])
    for orb1, orb2 in itertools.product(list(range(p.norb)), repeat=2):
        plt.plot(delta.data[:, 1, orb1, orb2].real)
        plt.plot(delta.data[:, 1, orb1, orb2].imag)

    plt.show()

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
            plot=False
            )

    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma

    test_symmetry_constraint_of_solve_eliashberg(g0_wk, gamma)
