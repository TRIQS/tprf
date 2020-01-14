# ----------------------------------------------------------------------

""" Test the functionality of the 'add_fake_bosonic_mesh' function """

# ----------------------------------------------------------------------

import numpy as np
import pytest

# ----------------------------------------------------------------------

from pytriqs.gf import Gf, MeshImFreq, MeshBrillouinZone, MeshProduct, Idx
from pytriqs.lattice import BrillouinZone, BravaisLattice

from triqs_tprf.lattice_utils import add_fake_bosonic_mesh

bz = BrillouinZone(BravaisLattice([[1,0],[0,1]]))
periodization_matrix = np.diag(np.array([10, 10, 1], dtype=np.int32))
bzmesh = MeshBrillouinZone(bz, periodization_matrix)

nmesh = MeshImFreq(beta=1, S='Fermion', n_max=1)

gf_nk = Gf(mesh=MeshProduct(nmesh, bzmesh), target_shape=(2,2))
gf_wnk = add_fake_bosonic_mesh(gf_nk)

np.testing.assert_allclose(gf_nk.data, gf_wnk[Idx(0), :, :].data)

gf_k = Gf(mesh=bzmesh, target_shape=(2,2))

with pytest.raises(ValueError):
    gf_wk = add_fake_bosonic_mesh(gf_k)

beta = 10
gf_wk = add_fake_bosonic_mesh(gf_k, beta=beta)

np.testing.assert_allclose(gf_k.data, gf_wk[Idx(0), :].data)
assert gf_wk.mesh[0].beta == beta
