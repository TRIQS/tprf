# ----------------------------------------------------------------------

import numpy as np

import triqs.utility.mpi as mpi

from triqs.gf import Gf, MeshReFreq, MeshBrZone
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice_utils import construct_densdens_V

# ----------------------------------------------------------------------

def test():
    nk = 10
    units = [(2,0,0), (0,2,0)]
    orb_pos = [(0,0,0), (0, 0.5, 0)]
    nb = len(orb_pos)

    bl = BravaisLattice(units=units, orbital_positions=orb_pos)
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, 1])
    V_k = construct_densdens_V(kmesh)

    # Check correct shape
    assert np.shape(V_k.data[:]) == (nk*nk,nb,nb,nb,nb)
    
    # Check for only density-density terms and symmetry
    for nb1 in range(nb):
        for nb2 in range(nb):
            for nb3 in range(nb):
                for nb4 in range(nb):
                    if ((nb1 == nb2) and (nb3 == nb4)):
                        assert np.allclose(V_k.data[:,nb1,nb2,nb3,nb4],
                              np.conjugate(V_k.data[:,nb3,nb4,nb1,nb2]))
                        continue

                    assert np.allclose(V_k.data[:,nb1,nb2,nb3,nb4], 0.0)

    # Check elements
    bl2 = BravaisLattice(units=units, orbital_positions=[(0,0,0)]*nb)
    bz2 = BrillouinZone(bl2)
    kmesh2 = MeshBrZone(bz2, [nk, nk, 1])
    V_k2 = construct_densdens_V(kmesh2)

    assert np.allclose(V_k2.data[:,0,0,0,0], V_k2.data[:,1,1,0,0])
    assert np.allclose(V_k2.data[:,0,0,0,0], V_k.data[:,0,0,0,0])
    assert np.allclose(V_k2.data[:,1,1,1,1], V_k.data[:,1,1,1,1])
    assert not np.allclose(V_k2.data[:,1,1,0,0], V_k.data[:,1,1,0,0])


if __name__ == "__main__":
    test()
