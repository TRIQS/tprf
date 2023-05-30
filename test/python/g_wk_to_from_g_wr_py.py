
import numpy as np

import triqs.utility.mpi as mpi

from triqs.gf import Gf, MeshImFreq, MeshProduct
from triqs.gf import MeshBrZone, MeshCycLat
from triqs.lattice import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_wk

bl = BravaisLattice([[1,0],[0,1]])
bz = BrillouinZone(bl)

lmesh = MeshCycLat(lattice=bl, dims=[10, 10, 1])
bzmesh = MeshBrZone(bz=bz, dims=[10, 10, 1])

e_k = Gf(mesh=bzmesh, target_shape=[1, 1])

for k in bzmesh:
    e_k[k] = -2*(np.cos(k[0]) + np.cos(k[1])) # does not work...
    
mesh = MeshImFreq(beta=1.0, S='Fermion', n_max=1024)
g0_wk = lattice_dyson_g0_wk(mu=1.0, e_k=e_k, mesh=mesh)

g0_wr = fourier_wk_to_wr(g0_wk)
g0_wk_ref = fourier_wr_to_wk(g0_wr)

np.testing.assert_array_almost_equal(g0_wk.data, g0_wk_ref.data)
