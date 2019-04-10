# ----------------------------------------------------------------------

""" Compare the naive implementation of the linearized Eliashberg product
and the one using Fourier transformations.
This test is quite computational intensive, because of the inefficiency of
the naive implementations. In the future the Fourier transformation 
implementation will subsitute the naive one.
"""

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, MeshImFreq, Idx

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk, solve_rpa_PH
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.lattice import eliashberg_product
from triqs_tprf.lattice import eliashberg_product_fft, eliashberg_product_fft_v2
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors

# ----------------------------------------------------------------------

def compare_deltas(deltas_1, deltas_2=None, static=False):

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

    i_max, j_max = diff.shape

    s = ""

    for i in range(i_max):

        for j in range(j_max):

            s += np.format_float_scientific(diff[i,j], precision=2, pad_left=3)
            s += "\t"

        s += "\n"
    print(s)






def compare_next_delta(p):

# -- Setup model, RPA susceptibilities and spin/charge interaction

    full_units = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    all_nn_hoppings = list(itertools.product([-1, 0, 1], repeat=p.dim)) 
    non_diagonal_hoppings = [ele for ele in all_nn_hoppings if sum(np.abs(ele)) == 1] 

    t = -p.t * np.eye(p.norbs)

    H = TBLattice(
                units = full_units[:p.dim],
                hopping = {hop : t for hop in non_diagonal_hoppings},
                orbital_positions = [(0,0,0)]*p.norbs,
                )

    e_k = H.on_mesh_brillouin_zone(n_k=[p.nk]*p.dim + [1]*(3-p.dim))

# A bigger w-mesh is needed to construct a Gamma with a twice as big w-mesh than GF

    wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw)
    wmesh_small = MeshImFreq(beta=p.beta, S='Fermion', n_max=int(p.small_factor*p.nw)+1)
    wmesh_big = MeshImFreq(beta=p.beta, S='Fermion', n_max=int(p.big_factor*p.nw)+1)

    g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)
    g0_wk_small = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh_small)
    g0_wk_big = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh_big)

    chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw)
    chi0_wk_small = imtime_bubble_chi0_wk(g0_wk_big, nw=int(p.small_factor*p.nw)+1)
    chi0_wk_big = imtime_bubble_chi0_wk(g0_wk_big, nw=int(p.big_factor*p.nw)+1)

    U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(p.norbs, p.U, 0, 0, 0)

    chi_s = solve_rpa_PH(chi0_wk, U_s)
    chi_c = solve_rpa_PH(chi0_wk, -U_c) # Minus for correct charge rpa equation
    chi_s_small = solve_rpa_PH(chi0_wk_small, U_s)
    chi_c_small = solve_rpa_PH(chi0_wk_small, -U_c) # Minus for correct charge rpa equation
    chi_s_big = solve_rpa_PH(chi0_wk_big, U_s)
    chi_c_big = solve_rpa_PH(chi0_wk_big, -U_c) # Minus for correct charge rpa equation

    gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)
    gamma_small = gamma_PP_singlet(chi_c_small, chi_s_small, U_c, U_s)
    gamma_big = gamma_PP_singlet(chi_c_big, chi_s_big, U_c, U_s)

# -- Preprocess gamma for the FFT implementations

    gamma_dyn_wk, gamma_const_k = split_into_dynamic_wk_and_constant_k(gamma)
    gamma_dyn_wk_small, gamma_const_k_small = split_into_dynamic_wk_and_constant_k(gamma_small)
    gamma_dyn_wk_big, gamma_const_k_big = split_into_dynamic_wk_and_constant_k(gamma_big)

    if not p.const:

        gamma.data[:] = gamma.data - p.U
        gamma_small.data[:] = gamma_small.data - p.U
        gamma_big.data[:] = gamma_big.data - p.U

        gamma_dyn_wk.data[:] = gamma.data
        gamma_dyn_wk_small.data[:] = gamma_small.data
        gamma_dyn_wk_big.data[:] = gamma_big.data

        gamma_const_k.data[:] = 0.0
        gamma_const_k_small.data[:] = 0.0
        gamma_const_k_big.data[:] = 0.0

    if not p.fit_const:

        gamma_dyn_wk.data[:] = gamma.data - p.U
        gamma_dyn_wk_small.data[:] = gamma_small.data - p.U
        gamma_dyn_wk_big.data[:] = gamma_big.data - p.U

        gamma_const_k.data[:] = p.U
        gamma_const_k_small.data[:] = p.U
        gamma_const_k_big.data[:] = p.U

    gamma_dyn_tr, gamma_const_r = dynamic_and_constant_to_tr(gamma_dyn_wk, gamma_const_k)
    gamma_dyn_tr_small, gamma_const_r_small = dynamic_and_constant_to_tr(gamma_dyn_wk_small,
                                                                            gamma_const_k_small)
    gamma_dyn_tr_big, gamma_const_r_big = dynamic_and_constant_to_tr(gamma_dyn_wk_big,
                                                                            gamma_const_k_big)

    # -- Creating Semi-Random input Delta

    np.random.seed(1337)

    v0 = g0_wk.copy()
    v0.data[:] = v0.data.real
    random_data = np.random.random(v0.data.shape[1:])
    freq_data = np.mean(np.abs(v0.data), axis=tuple(range(len(v0.data.shape))[1:]))
    not_randomized = 40
    start, stop = not_randomized, v0.data.shape[0]-not_randomized
    freq_data[start:stop] *= np.random.random(stop-start)

    v0.data[:] = np.tensordot(freq_data, random_data, axes=0)


    p.v0 = v0

# -- Test the Eliashberg equation
    
    print('summation')
    next_delta = eliashberg_product(gamma_big, g0_wk, p.v0)
    print('fft')
    next_delta_fft = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, p.v0)
    print('fft_small')
    next_delta_fft_small = eliashberg_product_fft(gamma_dyn_tr_small, gamma_const_r_small, g0_wk, p.v0)
    print('fft_big')
    next_delta_fft_big = eliashberg_product_fft(gamma_dyn_tr_big, gamma_const_r_big, g0_wk, p.v0)
    print('fft_v2')
    next_delta_fft_v2 = eliashberg_product_fft(gamma_dyn_tr, gamma_const_r, g0_wk, p.v0)


    from pytriqs.plot.mpl_interface import oplot, plt
    subp = [2, 6, 1]
    fig = plt.figure(figsize=(18, 15))

    deltas = [v0, next_delta, next_delta_fft, next_delta_fft_small, next_delta_fft_big, next_delta_fft_v2]
    titles = ['input', 'summation', 'fft', 'fft_small', 'fft_big', 'fft_v2']

    for k_point in [Idx(0,0,0), Idx(1,0,0)]:

        
        for delta, title in zip(deltas, titles):

            ax = plt.subplot(*subp); subp[-1] += 1
            oplot(delta[:, k_point])
            plt.title(title)

            ax.legend_ = None


    plt.show()

    diff = compare_deltas(deltas[1:])

    print_diff(diff)

    return deltas


    

# ----------------------------------------------------------------------

p = ParameterCollection(
        dim = 1,
        norbs = 1,
        t = 2.0,
        mu = 0.0,
        beta = 5,
        U = 1.0,
        nk = 4,
        nw = 200,
        const = True,
        fit_const = False,
        big_factor = 2,
        small_factor = 1.0,
        )

deltas_1 = compare_next_delta(p)

exit()

p.nw = 100


deltas_2 = compare_next_delta(p)

diff = compare_deltas(deltas_1[1:], deltas_2[1:], static=True)

print_diff(diff)


