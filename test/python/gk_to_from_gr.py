
import numpy as np

from pytriqs.gf import Gf, MeshImFreq, MeshProduct
from pytriqs.gf import MeshBrillouinZone, MeshCyclicLattice
from pytriqs.lattice import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import g0k_from_ek
from triqs_tprf.lattice import gr_from_gk
from triqs_tprf.lattice import gk_from_gr

bz = BrillouinZone(BravaisLattice([[1,0],[0,1]]))

periodization_matrix = np.diag(np.array([10, 10, 1], dtype=np.int32))

bzmesh = MeshBrillouinZone(bz, periodization_matrix)
lmesh = MeshCyclicLattice(bz.lattice, periodization_matrix)

ek = Gf(mesh=bzmesh, target_shape=[1, 1])

for idx, k in enumerate(bzmesh):
    ek[k] = -2*(np.cos(k[0]) + np.cos(k[1])) # does not work...
    #ek.data[idx] = -2*(np.cos(k[0]) + np.cos(k[1]))
    
mesh = MeshImFreq(beta=1.0, S='Fermion', n_max=1024)
g0k = g0k_from_ek(mu=1.0, ek=ek, mesh=mesh)

g0r = gr_from_gk(g0k, lmesh)
g0k_ref = gk_from_gr(g0r, bzmesh)

np.testing.assert_array_almost_equal(g0k.data, g0k_ref.data)
