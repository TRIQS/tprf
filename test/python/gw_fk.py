# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_fk
from triqs_tprf.lattice import lattice_dyson_g_fk
from triqs_tprf.lattice import lindhard_chi00_fk
from triqs_tprf.lattice import gw_sigma_fk_g0w0_spectral, gw_sigma_k_g0w0

from triqs.gf import Gf, MeshReFreq, inverse
from triqs.gf.mesh_product import MeshProduct

# ----------------------------------------------------------------------

nw = 3 #101
wmin = -5.0
wmax = 5.0
nk = 5
norb = 1
beta = 10.0
V = 1.0
mu = 0.0
delta = 0.001

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
fmesh = MeshReFreq(wmin, wmax, nw)
g_fk = lattice_dyson_g0_fk(mu=mu, e_k=e_k, mesh=fmesh, delta=delta)

V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
V_k.data[:] = V

print('--> pi_bubble')
PI_fk = lindhard_chi00_fk(e_k=e_k, mesh=fmesh, beta=beta, mu=mu, delta=delta)

print('--> screened_interaction_W')
Wr_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=V_k.target_shape)
for f in fmesh:
    fi = f.linear_index
    Wr_fk.data[fi,:] = V_k.data / (1.0 - PI_fk.data[fi,:] * V_k.data)


print('--> gw_self_energy')
sigma_fk = gw_sigma_fk_g0w0_spectral(mu=mu, beta=beta, e_k=e_k, mesh=fmesh, Wr_fk=Wr_fk, v_k=V_k, delta=delta)
sigma_k = gw_sigma_k_g0w0(mu=mu, beta=beta, e_k=e_k, v_k=V_k)

# TODO: Why does this return NaNs?
print(sigma_fk.data)


#f0_ind = len(fmesh)//2
#assert np.allclose( list(fmesh.values())[f0_ind], 0.0 )
#
#np.testing.assert_array_almost_equal(sigma_fk.data[f0_ind,:], sigma_k.data[:])

print('--> lattice_dyson_g_wk')
g_fk = lattice_dyson_g_fk(mu, e_k, sigma_fk, delta)
