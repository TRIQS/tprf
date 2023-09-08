# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.utility import mpi
from triqs.gf import Gf, MeshBrZone
from triqs.gf.mesh_product import MeshProduct
#from triqs_tprf.tight_binding import TBLattice
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import gaussian_dos

from triqs_tprf.lattice_utils import backfold_k

# ----------------------------------------------------------------------

def test_flat_dos():
    print("Single band")
    a0 = 2.0
    nk = 50
    sigma = 0.2
    mstar = 0.4 * 0.1314
    mu = 1.0

    print("--> setup dispersion")
    bl = BravaisLattice(units=[(a0,0,0), (0,a0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, 1])
    eps_k = Gf(mesh=kmesh, target_shape=[1])
    for k in kmesh:
        kfolded = backfold_k(k.value, kmesh)
        knorm = np.linalg.norm(kfolded)
        eps_k[k] = knorm**2.0 / (2.0 * mstar)

    print('--> get dos')
    dos = gaussian_dos(eps_k, mu, sigma)
    print(dos)

    print('--> compare with analytic result')
    A = a0*a0
    dos_ref = A*mstar / (2.0*np.pi)
    print(dos_ref)

    assert np.isclose(dos, dos_ref, rtol=1e-2)


def test_multiband():
    print("Multi band")
    a0 = 2.0
    nk = 50
    sigma = 0.2
    mstar1 = 0.4 * 0.1314
    mstar2 = 0.2 * 0.1314
    mu = 1.0
    shift = 0.2
    norb = 2

    print("--> setup dispersion")
    bl = BravaisLattice(units=[(a0,0,0), (0,a0,0)], orbital_positions=[(0,0,0)*norb])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, 1])

    eps_k = Gf(mesh=kmesh, target_shape=[norb])
    eps_k.data[:] = 0.0
    for k in kmesh:
        ki = k.data_index
        kfolded = backfold_k(k.value, kmesh)
        knorm = np.linalg.norm(kfolded)
        eps_k.data[ki,0] = knorm**2.0 / (2.0 * mstar1)
        eps_k.data[ki,1] = knorm**2.0 / (2.0 * mstar2) + shift

    print('--> get dos')
    dos = gaussian_dos(eps_k, mu, sigma)
    print(dos.real)

    print('--> compare with analytic result')
    A = a0*a0
    dos_ref1 = A*mstar1 / (2.0*np.pi)
    dos_ref2 = A*mstar2 / (2.0*np.pi)
    dos_ref = dos_ref1 + dos_ref2
    print(dos_ref)

    assert np.allclose(dos, dos_ref, rtol=1e-2)


if __name__ == "__main__":
    test_flat_dos()
    test_multiband()
