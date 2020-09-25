# ----------------------------------------------------------------------

""" Compare the summation implementation of the linearized Eliashberg product
and the one using Fourier transformations.

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """

# ----------------------------------------------------------------------

import numpy as np
from triqs.plot.mpl_interface import oplot, plt
import warnings 

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs.gf import Idx
from triqs_tprf.utilities import create_eliashberg_ingredients
from triqs_tprf.lattice import eliashberg_product, eliashberg_product_fft
from triqs_tprf.eliashberg import semi_random_initial_delta, preprocess_gamma_for_fft

# ----------------------------------------------------------------------

def compare_deltas(deltas_1, deltas_2=None, static=False):
    """ Build comparison matrix of list of Gf
    """
    if not deltas_2:
        deltas_2 = deltas_1

    if static:
        deltas_1 = [ele[Idx(0), :] for ele in deltas_1]
        deltas_2 = [ele[Idx(0), :] for ele in deltas_2]

    diff = np.zeros(shape=(len(deltas_1), len(deltas_2)))

    for i, delta_1 in enumerate(deltas_1):
        for j, delta_2 in enumerate(deltas_2):

            diff[i,j] = np.max(np.abs(delta_1.data - delta_2.data))

    return diff

def print_diff(diff):
    """ Print output of 'compare_deltas' more readable
    """
    i_max, j_max = diff.shape

    s = ""
    s += "\n Differences matrix\n"
    dashes = "-"*14*diff.shape[0] + "\n"
    s += dashes

    for i in range(i_max):
        for j in range(j_max):
            try:
                s += np.format_float_scientific(diff[i,j], precision=2, pad_left=3)
            except AttributeError:
                s += str(diff[i,j])
            s += "\t"
        s += "\n"
    s += dashes
    print(s)

def test_eliashberg_product_for_same_initital_delta(g0_wk, gamma, gamma_big):
    initial_delta = semi_random_initial_delta(g0_wk, seed=1337)

    next_delta_summation = eliashberg_product(gamma_big, g0_wk, initial_delta)

    gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma)
    next_delta_fft = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, initial_delta)

    diff = compare_deltas([next_delta_summation, next_delta_fft])

    print_diff(diff)
    np.testing.assert_allclose(diff, 0, atol=p.atol)
    print('The summation and FFT implementation of the eliashberg product'
                                                ' both yield the same result.')
def test_eliashberg_product_for_different_initital_delta(g0_wk, gamma, gamma_big):
    initial_delta = semi_random_initial_delta(g0_wk, seed=1337)

    next_delta_summation = eliashberg_product(gamma_big, g0_wk, initial_delta)

    gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma)
    initial_delta = semi_random_initial_delta(g0_wk, seed=1338)
    next_delta_fft = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, initial_delta)

    diff = compare_deltas([next_delta_summation, next_delta_fft])

    print_diff(diff)
    try:
        np.testing.assert_allclose(diff, 0, atol=p.atol)
        raise ValueError
    except AssertionError:
        print('The summation and FFT implementation of the eliashberg product'
        ' both yield DIFFERENT results, as expected when using a different inital delta.')

def plot_output(g0_wk, gamma):
    initial_delta = semi_random_initial_delta(g0_wk)

    next_delta_summation = eliashberg_product(gamma_big, g0_wk, initial_delta)

    gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma)
    next_delta_fft = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, initial_delta)

    deltas = [initial_delta, next_delta_summation, next_delta_fft]

    warnings.filterwarnings("ignore") #ignore some matplotlib warnings
    subp = [4, 3, 1]
    fig = plt.figure(figsize=(18, 15))

    titles = ['Input', 'Summation', 'FFT']

    for k_point in [Idx(0,0,0), Idx(1,0,0)]:

        ax = plt.subplot(*subp); subp[-1] += 1
        oplot(g0_wk[:, k_point])
        plt.title('GF')

        ax = plt.subplot(*subp); subp[-1] += 1
        oplot(gamma[:, k_point])
        plt.title('Gamma')

        ax = plt.subplot(*subp); subp[-1] += 1
        oplot(gamma_dyn_tr[:, k_point])
        plt.title('Gamma dyn tr')

        for delta, title in zip(deltas, titles):

            ax = plt.subplot(*subp); subp[-1] += 1
            oplot(delta[:, k_point])
            plt.title(title)

            ax.legend_ = None

    plt.show()
    
#================================================================================ 

if __name__ == '__main__':
    p = ParameterCollection(
            dim = 1,
            norb = 2,
            t1 = 1.0,
            t2 = 0.5,
            t12 = 0.1,
            t21 = 0.1,
            mu = 0.0,
            beta = 1,
            U = 1.0,
            Up = 0.8,
            J = 0.1,
            Jp = 0.1,
            nk = 3,
            nw = 30,
            atol = 1e-8,
            )

    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma
    # For the eliashberg SUM procedure a Gamma with a twice as big w-mesh then the GF is needed.
    big_nw = 2*p.nw + 1
    eliashberg_ingredients_big = create_eliashberg_ingredients(p.alter(nw=big_nw))
    gamma_big = eliashberg_ingredients_big.gamma

    test_eliashberg_product_for_same_initital_delta(g0_wk, gamma, gamma_big)
    test_eliashberg_product_for_different_initital_delta(g0_wk, gamma, gamma_big)
    #plot_output(g0_wk, gamma)
