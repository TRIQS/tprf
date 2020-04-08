# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.gw import bubble_PI_wk
from triqs_tprf.gw import retarded_screened_interaction_Wr_wk
from triqs_tprf.gw import gw_sigma_wk

from pytriqs.gf import Gf, MeshImFreq, Idx

# ----------------------------------------------------------------------

nw = 100
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

e_k = t_r.on_mesh_brillouin_zone(n_k=(8, 1, 1))

print(e_k.data)

kmesh = e_k.mesh
wmesh = MeshImFreq(beta, 'Fermion', nw)
g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
V_k.data[:] = V

print('--> pi_bubble')
PI_wk = bubble_PI_wk(g_wk)

print('--> screened_interaction_W')
Wr_wk = retarded_screened_interaction_Wr_wk(PI_wk, V_k)

print('--> gw_self_energy')
sigma_wk = gw_sigma_wk(Wr_wk, g_wk, fft_flag=True)
sigma_wk_ref = gw_sigma_wk(Wr_wk, g_wk, fft_flag=False)
np.testing.assert_array_almost_equal(sigma_wk.data, sigma_wk_ref.data)

print('--> lattice_dyson_g_wk')
g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
