# ----------------------------------------------------------------------

import itertools
from itertools import product
from tqdm import tqdm
import numpy as np

# ----------------------------------------------------------------------

from triqs.utility import mpi
from triqs.gf import Gf, MeshBrZone
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import densdens_V_orb_to_band_basis

# ----------------------------------------------------------------------

def orb_to_band_python(V_orb_q, psi_k):
    kmesh = psi_k.mesh
    nb = np.shape(V_orb_q.data[:])[-1]
    V_band_kkp = Gf(mesh=MeshProduct(kmesh, kmesh), target_shape=[nb, nb])
    V_band_kkp.data[:] = 0.0

    for k in kmesh:
        ki = k.data_index
        for kp in kmesh:
            kpi = kp.data_index
            V_orb_kmkp = V_orb_q(k-kp)

            for a,b in product(range(nb), range(nb)):
                for i,j in product(range(nb), range(nb)):
                    V_band_kkp.data[ki,kpi,i,j] += \
                    np.conj(psi_k.data[ki,a,i]) * \
                    np.conj(psi_k.data[kpi,b,j]) *\
                    psi_k.data[ki,b,j] *\
                    psi_k.data[kpi,a,i] * V_orb_kmkp[a,a,b,b]

    return V_band_kkp


def test_orb_to_band_basis():
    print("Test against python implementation")
    nk = 5
    norb = 3

    print("--> setup kmesh")
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, 1])

    print("--> setup V in orbital basis")
    V_orb_q = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_orb_q.data[:] = np.random.rand(*np.shape(V_orb_q.data[:]))

    print('--> setup psi in orbital basis')
    psi_k = Gf(mesh=kmesh, target_shape=[norb]*2)
    psi_k.data[:] = np.random.rand(*np.shape(psi_k.data[:]))

    print('--> get V in band basis (c++)')
    V_band_kkp = densdens_V_orb_to_band_basis(V_orb_q, psi_k)

    print('--> get V in band basis (python)')
    V_band_kkp_ref = orb_to_band_python(V_orb_q, psi_k)

    print('--> compare')
    assert np.allclose(V_band_kkp.data[:], V_band_kkp_ref.data[:])

def test_trivial():
    print("Trivial test")
    nk = 1
    norb = 1

    print("--> setup kmesh")
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, 1])

    print("--> setup V in orbital basis")
    V_orb_q = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_orb_q.data[:] = np.random.rand(*np.shape(V_orb_q.data[:]))

    print('--> setup psi in orbital basis')
    psi_k = Gf(mesh=kmesh, target_shape=[norb]*2)
    psi_k.data[:] = np.random.rand(*np.shape(psi_k.data[:]))

    print('--> get V in band basis (c++)')
    V_band_kkp = densdens_V_orb_to_band_basis(V_orb_q, psi_k)

    print('--> compare against trivial solution')
    V_band_kkp_ref = Gf(mesh=MeshProduct(kmesh,kmesh), target_shape=[norb]*2)
    V_band_kkp_ref.data[:,0,0] = psi_k.data[:,0,0]**2.0 * np.conj(psi_k.data[:,0,0])**2.0 * V_orb_q.data[:,0,0,0,0]
    assert np.allclose(V_band_kkp.data[:], V_band_kkp_ref.data[:])


if __name__ == "__main__":
    test_trivial()
    test_orb_to_band_basis()
