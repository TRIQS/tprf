# ----------------------------------------------------------------------

""" Test the functionality of the 'add_fake_bosonic_mesh' function """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import Gf, MeshImFreq, MeshBrillouinZone, MeshProduct, Idx
from triqs.lattice import BrillouinZone, BravaisLattice

from triqs_tprf.lattice_utils import add_fake_bosonic_mesh


def test_add_fake_bosonic_mesh_with_gf_nk(bzmesh):
    nmesh = MeshImFreq(beta=1, S="Fermion", n_max=1)

    gf_nk = Gf(mesh=MeshProduct(nmesh, bzmesh), target_shape=(2, 2))
    gf_wnk = add_fake_bosonic_mesh(gf_nk)

    np.testing.assert_allclose(gf_nk.data, gf_wnk[Idx(0), :, :].data)


def test_add_fake_bosonic_mesh_with_gf_k_without_beta(bzmesh):
    gf_k = Gf(mesh=bzmesh, target_shape=(2, 2))

    try:
        gf_wk = add_fake_bosonic_mesh(gf_k)
    except ValueError:
        pass


def test_add_fake_bosonic_mesh_with_gf_k_with_beta(bzmesh):
    gf_k = Gf(mesh=bzmesh, target_shape=(2, 2))
    beta = 10
    gf_wk = add_fake_bosonic_mesh(gf_k, beta=beta)

    np.testing.assert_allclose(gf_k.data, gf_wk[Idx(0), :].data)
    assert gf_wk.mesh[0].beta == beta


if __name__ == "__main__":
    bz = BrillouinZone(BravaisLattice([[1, 0], [0, 1]]))
    periodization_matrix = np.diag(np.array([10, 10, 1], dtype=np.int32))
    bzmesh = MeshBrillouinZone(bz, periodization_matrix)

    test_add_fake_bosonic_mesh_with_gf_nk(bzmesh)
    test_add_fake_bosonic_mesh_with_gf_k_without_beta(bzmesh)
    test_add_fake_bosonic_mesh_with_gf_k_with_beta(bzmesh)
