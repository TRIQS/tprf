# ----------------------------------------------------------------------

""" Test the symmetrizing feature of the solve_eliashberg function
"""
# ----------------------------------------------------------------------

import itertools
import functools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq, Idx

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors
from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------

import matplotlib.pyplot as plt

# ----------------------------------------------------------------------

from triqs_tprf.symmetries import enforce_symmetry, check_symmetry

# ----------------------------------------------------------------------

p = ParameterCollection(
        dim = 2,
        norbs = 2,
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

# -- Setup model, RPA susceptibilities, spin/charge interaction and gamma
full_units = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
all_nn_hoppings = list(itertools.product([-1, 0, 1], repeat=p.dim)) 
non_diagonal_hoppings = [ele for ele in all_nn_hoppings if sum(np.abs(ele)) == 1] 

# -- Create hopping matrix for two-band model
t = -np.array([[p.t1, p.t12], [p.t21, p.t2]])

H = TBLattice(
            units = full_units[:p.dim],
            hopping = {hop : t for hop in non_diagonal_hoppings},
            orbital_positions = [(0,0,0)]*p.norbs,
            )

e_k = H.on_mesh_brillouin_zone(n_k=[p.nk]*p.dim + [1]*(3-p.dim))

wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw)
g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)

chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw)

U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(p.norbs, p.U, p.Up,
                                                                        p.J, p.Jp)

chi_s = solve_rpa_PH(chi0_wk, U_s)
chi_c = solve_rpa_PH(chi0_wk, -U_c) # Minus for correct charge rpa equation

gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)

# -- Test symmetrizing function on eliashberg
variables=["frequency", "momentum", "orbital"]

# Use all combinations
symmetry_set = list(itertools.product(["even", "odd"], repeat=3))

translate_symmetries = {"even" : +1, "odd" : -1, None : None}

for symmetries in symmetry_set:
    symmetrize = functools.partial(enforce_symmetry, 
                                   variables=variables,
                                   symmetries=symmetries)

    E, eigen_modes = solve_eliashberg(gamma, g0_wk, product='FFT', solver='IRAM',
                                        symmetrize_fct=symmetrize)

    expected_symmetries = {variable : translate_symmetries[symmetry] \
                    for (variable, symmetry) in zip(variables, symmetries)}

    for delta in eigen_modes:
        produced_symmetries = check_symmetry(delta)
        if not expected_symmetries == produced_symmetries:
            raise AssertionError("Incorrect symmetries were produced.")

        if p.plot:
            fig, axes = plt.subplots(3, 3)

            vmax = np.max(np.abs(delta[Idx(0),:].data))

            for orb1, orb2 in itertools.product(range(p.norbs), repeat=2):
                shape = (p.nk, p.nk, p.norbs, p.norbs)
                data = delta[Idx(0), :].data.reshape(shape)
                plt.sca(axes[orb1,orb2])
                plt.imshow(data[:,:,orb1,orb2].real, cmap="RdBu_r",
                                                         vmax=vmax, vmin=-vmax)
                plt.colorbar()

            plt.sca(axes[-1,-1])
            for orb1, orb2 in itertools.product(range(p.norbs), repeat=2):
                plt.plot(delta.data[:, 10, orb1, orb2].real)
                plt.plot(delta.data[:, 10, orb1, orb2].imag)

            plt.show()

