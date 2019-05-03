# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr
from triqs_tprf.lattice import chi0_tr_from_grt_PH
from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_wk_from_chi_wr

from triqs_tprf.lattice import screened_interaction_W
from triqs_tprf.lattice import gw_self_energy

from pytriqs.gf import Gf, MeshImFreq, Idx

# ----------------------------------------------------------------------
def pi_bubble(g_wk):

    nw = len(g_wk.mesh.components[0]) / 2
    g_wr = fourier_wk_to_wr(g_wk)
    g_tr = fourier_wr_to_tr(g_wr)
    del g_wr

    PI_tr = chi0_tr_from_grt_PH(g_tr)
    del g_tr
    PI_wr = chi_wr_from_chi_tr(PI_tr, nw=nw)
    del PI_tr
    PI_wk = chi_wk_from_chi_wr(PI_wr)
    del PI_wr

    return PI_wk    
    
# ----------------------------------------------------------------------

nw = 40
norb = 1
beta = 10.0
V = 10.0
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

e_k = t_r.on_mesh_brillouin_zone(n_k=(1, 1, 1))

kmesh = e_k.mesh
wmesh = MeshImFreq(beta, 'Fermion', nw)
g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
V_k.data[:] = V

print('--> pi_bubble')
PI_wk = pi_bubble(g_wk)

print('--> screened_interaction_W')
W_wk = screened_interaction_W(PI_wk, V_k)

diff = np.max(np.abs(W_wk.data - V))
print 'W - V =', diff
np.testing.assert_array_almost_equal(W_wk.data, V + np.zeros_like(W_wk.data), decimal=2)

print('--> gw_self_energy')
sigma_wk = gw_self_energy(W_wk, g_wk)

print('--> lattice_dyson_g_wk')
g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
