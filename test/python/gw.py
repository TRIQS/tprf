# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.gw import bubble_PI_wk
from triqs_tprf.gw import dynamical_screened_interaction_W
from triqs_tprf.gw import gw_sigma
from triqs_tprf.gw import g0w_sigma

from triqs.gf import Gf, MeshImFreq, Idx

# ----------------------------------------------------------------------

nw = 1000
nk = 8
norb = 1
beta = 10.0
V = 1.0
mu = 0.0

t = -1.0 * np.eye(norb)

t_r = TBLattice(
    units = [(1, 0, 0)],
    hopping = {
        (+1,) : t,
        (-1,) : t,
        },
    orbital_positions = [(0,0,0)]*norb,
    )

kmesh = t_r.get_kmesh(n_k=(nk, 1, 1))
e_k = t_r.fourier(kmesh)

kmesh = e_k.mesh
wmesh = MeshImFreq(beta, 'Fermion', nw)
g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
V_k.data[:] = V

print('--> pi_bubble')
PI_wk = bubble_PI_wk(g0_wk)

print('--> screened_interaction_W (static bare interaction)')
Wr_wk = dynamical_screened_interaction_W(PI_wk, V_k)

Wr_full_wk = Gf(mesh=Wr_wk.mesh, target_shape=[norb]*4)
for w in Wr_full_wk.mesh.components[0]:
    Wr_full_wk[w,:] = Wr_wk[w,:] + V_k

print('--> gw_self_energy')
sigma_wk = gw_sigma(Wr_full_wk, g0_wk)

sigma_dyn_wk = gw_sigma(Wr_wk, g0_wk)
sigma_stat_k = gw_sigma(V_k, g0_wk)
sigma_wk_ref = Gf(mesh=sigma_dyn_wk.mesh, target_shape=sigma_dyn_wk.target_shape)
for w in wmesh:
    sigma_wk_ref[w,:] = sigma_dyn_wk[w,:] + sigma_stat_k

np.testing.assert_array_almost_equal(sigma_wk.data, sigma_wk_ref.data, decimal=5)





sigma_k = gw_sigma(V_k, g0_wk)
sigma_k_ref = g0w_sigma(mu=mu, beta=beta, e_k=e_k, v_k=V_k)

np.testing.assert_array_almost_equal(sigma_k.data, sigma_k_ref.data, decimal=4)




print('--> lattice_dyson_g_wk')
g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
