# ----------------------------------------------------------------------

""" Test if memory saving version of the bare bubble chi0 is implemented
correctly. """

# ----------------------------------------------------------------------

import itertools
import time

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------

p = ParameterCollection(
        dim = 2,
        norbs = 2,
        t = 1.0,
        mu = 0.0,
        beta = 1,
        nk = 16,
        nw_g = 100,
        nw_chi = 10,
        )

# -- Setup model
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

wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw_g)

g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)

# -- Calculate bare bubble 
chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw_chi, save_memory=False)
chi0_wk_save_memory = imtime_bubble_chi0_wk(g0_wk, nw=p.nw_chi, save_memory=True)

assert np.allclose(chi0_wk.data, chi0_wk_save_memory.data)
