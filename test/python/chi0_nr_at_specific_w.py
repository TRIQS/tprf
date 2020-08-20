# ----------------------------------------------------------------------

""" Test if accessing chi0_wnr only a a specififc w point gives the same
result as the function used for a whole bosonic Matsubara mesh. """


# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq, Idx, MeshProduct

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi0r_from_gr_PH, chi0_nr_from_gr_PH_at_specific_w

# ----------------------------------------------------------------------

p = ParameterCollection(
        dim = 2,
        norbs = 2,
        t = 1.0,
        mu = 0.0,
        beta = 1,
        nk = 16,
        nw_g = 100,
        nw_chi = 20,
        nwf = 10,
        nw_index = 3,
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
g0_wr = fourier_wk_to_wr(g0_wk)

chi0_wnr = chi0r_from_gr_PH(nw=p.nw_chi, nn=p.nwf, g_nr=g0_wr)
chi0_nr_at_specific_w = chi0_nr_from_gr_PH_at_specific_w(nw_index=p.nw_index, nn=p.nwf, g_nr=g0_wr)

assert isinstance(chi0_nr_at_specific_w.mesh, MeshProduct)
assert isinstance(chi0_nr_at_specific_w.mesh[0], MeshImFreq)
np.testing.assert_allclose(chi0_wnr[Idx(p.nw_index), :, :].data, chi0_nr_at_specific_w.data)


