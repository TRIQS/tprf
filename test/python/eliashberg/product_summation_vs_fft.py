# ----------------------------------------------------------------------

""" Compare the summation implementation of the linearized Eliashberg product
and the one using Fourier transformations.
"""

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, MeshImFreq, Idx
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

            s += np.format_float_scientific(diff[i,j], precision=2, pad_left=3)
            s += "\t"

        s += "\n"
    s += dashes
    print(s)

def compare_next_delta(p):
    eliashberg_ingredients = create_eliashberg_ingredients(p)
    g0_wk = eliashberg_ingredients.g0_wk
    gamma = eliashberg_ingredients.gamma
    U_c = eliashberg_ingredients.U_c
    U_s = eliashberg_ingredients.U_s

    ## A bigger w-mesh is needed to construct a Gamma with a twice as big w-mesh than GF
    big_nw = 2*p.nw + 1
    eliashberg_ingredients_big = create_eliashberg_ingredients(p.alter(nw=big_nw))
    gamma_big = eliashberg_ingredients_big.gamma

    # -- Preprocess gamma for the FFT implementations
    if p.fit_const:
        gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma, None)
    else:
        gamma_dyn_tr, gamma_const_r = preprocess_gamma_for_fft(gamma, 0.5*(U_s + U_c))

    # -- Creating Semi-Random input Delta
    v0 = semi_random_initial_delta(g0_wk, nr_factor=p.nr_factor, seed=1337)
    p.v0 = v0

    # -- Test the Eliashberg product
    print('Start the summation')
    next_delta = eliashberg_product(gamma_big, g0_wk, p.v0)
    print('Start the FFT')
    next_delta_fft = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, p.v0)

    deltas = [v0, next_delta, next_delta_fft]

    if p.plot:

        from pytriqs.plot.mpl_interface import oplot, plt
        import warnings 
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

    diff = compare_deltas(deltas[1:])

    print_diff(diff)
    try:
        np.testing.assert_allclose(diff, 0, atol=p.atol)
    except AssertionError as e:
        print('The test failed for the parameter set:')
        p.__dict__.pop("v0")
        print(p)
        raise e

    return deltas
    
#================================================================================ 

if __name__ == '__main__':

    p = ParameterCollection(
            dim = 1,
            norb = 1,
            t = 2.0,
            mu = 0.0,
            beta = 5,
            U = 1.0,
            Up = 0.8,
            J = 0.1,
            Jp = 0.1,
            nk = 3,
            nw = 150,
            nr_factor = 0.5,
            fit_const = False,
            big_factor = 2,
            atol = 1e-8,
            plot = False,
            )

    for norb in [1, 2]:
        p.norb = norb
        deltas = compare_next_delta(p)

    print('The summation and FFT implementation of the eliashberg product'
                                                ' both yield the same result.')

    deltas_with_fit = compare_next_delta(p.alter(fit_const=True))

    diff = compare_deltas(deltas[2:], deltas_with_fit[2:])

    print('Compare explicit given constant vs. fit:')

    print_diff(diff)
    np.testing.assert_allclose(diff, 0, atol=p.atol)

    print('Fitting the constant part works.')
