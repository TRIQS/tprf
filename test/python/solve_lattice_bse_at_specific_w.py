# ----------------------------------------------------------------------

""" Test if calculating the lattice susceptibility via the Bethe-
    Salpeter equation for a specific \omega gives the same result as the
    function used for a whole bosonic Matsubara mesh. """


# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import Gf, MeshImFreq, Idx, MeshProduct, MeshBrillouinZone

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.bse import solve_lattice_bse, solve_lattice_bse_at_specific_w

# ----------------------------------------------------------------------

p = ParameterCollection(
    dim=2,
    norbs=2,
    t=1.0,
    mu=0.0,
    beta=1,
    nk=16,
    nw_g=100,
    nw_gamma=10,
    nwf=10,
    nw_index=3,
    seed=101,
)

# -- Setup model
full_units = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
all_nn_hoppings = list(itertools.product([-1, 0, 1], repeat=p.dim))
non_diagonal_hoppings = [ele for ele in all_nn_hoppings if sum(np.abs(ele)) == 1]

t = -p.t * np.eye(p.norbs)

H = TBLattice(
    units=full_units[: p.dim],
    hopping={hop: t for hop in non_diagonal_hoppings},
    orbital_positions=[(0, 0, 0)] * p.norbs,
)

e_k = H.on_mesh_brillouin_zone(n_k=[p.nk] * p.dim + [1] * (3 - p.dim))

wmesh = MeshImFreq(beta=p.beta, S="Fermion", n_max=p.nw_g)

g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)

wmesh_gamma = MeshImFreq(beta=p.beta, S="Boson", n_max=p.nw_gamma)
nmesh_gamma = MeshImFreq(beta=p.beta, S="Fermion", n_max=p.nwf)

gamma_wnn = Gf(
    mesh=MeshProduct(wmesh_gamma, nmesh_gamma, nmesh_gamma),
    target_shape=2 * g0_wk.target_shape,
)

np.random.seed(p.seed)
gamma_wnn.data[:] = np.random.rand(*gamma_wnn.data.shape)

chi_kw, chi0_kw = solve_lattice_bse(g0_wk, gamma_wnn)
chi_k_at_specific_w, chi0_k_at_specific_w = solve_lattice_bse_at_specific_w(
    g0_wk, gamma_wnn, nw_index=p.nw_index
)

assert isinstance(chi_k_at_specific_w.mesh, MeshBrillouinZone)

np.testing.assert_allclose(chi0_kw[:, Idx(p.nw_index)].data, chi0_k_at_specific_w.data, atol=10e-16)
np.testing.assert_allclose(chi_kw[:, Idx(p.nw_index)].data, chi_k_at_specific_w.data, atol=10e-16)
