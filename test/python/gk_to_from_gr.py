
import numpy as np

from pytriqs.gf import Gf, MeshImFreq, MeshProduct, MeshBrillouinZone
from pytriqs.lattice import BrillouinZone, BravaisLattice

from pytriqs.applications.tprf.lattice import g0k_from_ek
from pytriqs.applications.tprf.lattice import gr_from_gk
from pytriqs.applications.tprf.lattice import gk_from_gr

bz = BrillouinZone(BravaisLattice([[1,0],[0,1]]))
bzmesh = MeshBrillouinZone(bz, n_k=10)
ek = Gf(mesh=bzmesh, target_shape=[1, 1])

for idx, k in enumerate(bzmesh):
    #ek[k] = -2*t*(np.cos(k[0]) + np.cos(k[1])) # does not work...
    ek.data[idx] = -2*(np.cos(k[0]) + np.cos(k[1]))

mesh = MeshImFreq(beta=1.0, S='Fermion', n_max=1024)
g0k = g0k_from_ek(mu=1.0, ek=ek, mesh=mesh)

g0r = gr_from_gk(g0k)
g0k_ref = gk_from_gr(g0r, bz)

np.testing.assert_array_almost_equal(g0k.data, g0k_ref.data)
