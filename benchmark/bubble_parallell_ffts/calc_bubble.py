# ----------------------------------------------------------------------

import time
import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq

from triqs_tprf.ParameterCollection import ParameterCollection

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------

dim = 3
norbs = 6
nk = 8
nw = 128
beta = 40.0

full_units = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
all_nn_hoppings = list(itertools.product([-1, 0, 1], repeat=dim)) 
non_diagonal_hoppings = [ele for ele in all_nn_hoppings if sum(np.abs(ele)) == 1] 

t = -1.0 * np.eye(norbs)

H = TBLattice(
            units = full_units[:dim],
            hopping = {hop : t for hop in non_diagonal_hoppings},
            orbital_positions = [(0,0,0)]*norbs,
            )
e_k = H.on_mesh_brillouin_zone(n_k=[nk]*dim + [1]*(3-dim))

wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)

print('--> Dyson g0_wk')
t = time.time()
g0_wk = lattice_dyson_g0_wk(mu=0.0, e_k=e_k, mesh=wmesh)
print(' {} seconds'.format(time.time() - t))

print('--> Bubble chi0_wk')
t = time.time()
chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=nw)
print(' {} seconds'.format(time.time() - t))

