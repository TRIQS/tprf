import numpy as np
import sys
import time
import os

import matplotlib.pyplot as plt

from triqs.gfs import *
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk, dynamic_and_constant_to_tr

from triqs_tprf.lattice import eliashberg_F_wk, eliashberg_g_wk, eliashberg_gbar_wk

def test_multiorbital_F_wk():
    beta = 10.0
    nw = 100
    nk = 1
    decimal = 10

    h = np.array([[0.4, 0.3],
                  [0.3,-0.7]])
    hbar = np.array([[-0.2, 0.1],
                     [ 0.1, 0.3]])
    delta = np.array([[0.2,  0.05],
                      [0.05, 0.5 ]])

    print('--> Construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, 1, 1])

    wmesh_f  = MeshImFreq(beta, 'Fermion', nw)

    print('--> Construct BdG Hamiltonian')
    H_BdG = np.block([
        [h,              delta],
        [delta.conj().T, hbar ]
    ])

    assert np.allclose(H_BdG, np.conjugate(np.transpose(H_BdG)))

    Norb = np.shape(h)[0]
    G_BdG_wk = Gf(mesh=MeshProduct(wmesh_f, kmesh), target_shape=[2*Norb]*2)
    for k in kmesh:
        G_BdG_wk[:,k] << inverse(iOmega_n - H_BdG)

    g_BdG_wk = G_BdG_wk[:Norb,:Norb]
    gbar_BdG_wk = G_BdG_wk[Norb:,Norb:]
    F_BdG_wk = -G_BdG_wk[:Norb,Norb:] # Negative sign due to convention of TPRF
    Fstar_BdG_wk = -G_BdG_wk[Norb:,:Norb] 

    print("--> Construct Greens functions")
    h_k = Gf(mesh=kmesh, target_shape=[Norb]*2)
    h_k.data[:] = h
    hbar_k = Gf(mesh=kmesh, target_shape=[Norb]*2)
    hbar_k.data[:] = hbar

    g0_wk = lattice_dyson_g0_wk(mu=0.0, e_k=h_k, mesh=wmesh_f)
    gbar0_wk = lattice_dyson_g0_wk(mu=0.0, e_k=hbar_k, mesh=wmesh_f)

    delta_wk = Gf(mesh=MeshProduct(wmesh_f, kmesh), target_shape=[Norb]*2)
    for w in wmesh_f:
        for k in kmesh:
            delta_wk[w,k] = delta

    print('--> Compare F_wk to BdG')
    F_wk = eliashberg_F_wk(g0_wk, gbar0_wk, delta_wk)

    print("F:", F_wk.data[nw,0,:,:])
    print("F_BdG:", F_BdG_wk.data[nw,0,:,:])
    print("Fstar_BdG:", Fstar_BdG_wk.data[nw,0,:,:])

    np.testing.assert_array_almost_equal(F_wk.data[:], 
                                         F_BdG_wk.data[:], 
                                         decimal=decimal)

    print('--> Compare electron g_wk to BdG')
    g_wk = eliashberg_g_wk(g0_wk, gbar0_wk, delta_wk)
   
    print("g:", g_wk.data[nw,0,:,:])
    print("g_BdG:", g_BdG_wk.data[nw,0,:,:])

    np.testing.assert_array_almost_equal(g_wk.data[:], 
                                         g_BdG_wk.data[:], 
                                         decimal=decimal)

    print('--> Compare hole gbar_wk to BdG')
    gbar_wk = eliashberg_gbar_wk(g0_wk, gbar0_wk, delta_wk)
   
    print("gbar:", gbar_wk.data[nw,0,:,:])
    print("gbar_BdG:", gbar_BdG_wk.data[nw,0,:,:])

    np.testing.assert_array_almost_equal(gbar_wk.data[:], 
                                         gbar_BdG_wk.data[:], 
                                         decimal=decimal)

    print('--> Construct BdG Hamiltonian when hbar=-h.T')
    H_BdG = np.block([
        [h,              delta],
        [delta.conj().T, -h.T ]
    ])

    assert np.allclose(H_BdG, np.conjugate(np.transpose(H_BdG)))

    Norb = np.shape(h)[0]
    G_BdG_wk = Gf(mesh=MeshProduct(wmesh_f, kmesh), target_shape=[2*Norb]*2)
    for k in kmesh:
        G_BdG_wk[:,k] << inverse(iOmega_n - H_BdG)

    g_BdG_wk = G_BdG_wk[:Norb,:Norb]
    gbar_BdG_wk = G_BdG_wk[Norb:,Norb:]
    F_BdG_wk = -G_BdG_wk[:Norb,Norb:] # Negative sign due to convention of TPRF
    Fstar_BdG_wk = -G_BdG_wk[Norb:,:Norb] 

    print('--> Compare to BdG when gbar_wk is not supplied')
    g0_wk = lattice_dyson_g0_wk(mu=0.0, e_k=h_k, mesh=wmesh_f)

    F_wk = eliashberg_F_wk(g0_wk, delta_wk)
    g_wk = eliashberg_g_wk(g0_wk, delta_wk)
    gbar_wk = eliashberg_gbar_wk(g0_wk, delta_wk)

    print('--> F_wk')
    print("F:", F_wk.data[nw,0,:,:])
    print("F_BdG:", F_BdG_wk.data[nw,0,:,:])
    print("Fstar_BdG:", Fstar_BdG_wk.data[nw,0,:,:])

    np.testing.assert_array_almost_equal(F_wk.data[:], 
                                         F_BdG_wk.data[:], 
                                         decimal=decimal)

    print('--> g_wk')
    print("g:", g_wk.data[nw,0,:,:])
    print("g_BdG:", g_BdG_wk.data[nw,0,:,:])

    np.testing.assert_array_almost_equal(g_wk.data[:], 
                                         g_BdG_wk.data[:], 
                                         decimal=decimal)

    print('--> gbar_wk')
    print("gbar:", gbar_wk.data[nw,0,:,:])
    print("gbar_BdG:", gbar_BdG_wk.data[nw,0,:,:])

    np.testing.assert_array_almost_equal(gbar_wk.data[:], 
                                         gbar_BdG_wk.data[:], 
                                         decimal=decimal)

    print('--> Check G(w) = -Gbar(-w)')
    g_ref = -gbar_wk.conjugate().data[:]
    g_ref = np.swapaxes(g_ref, 2, 3)
    np.testing.assert_array_almost_equal(g_wk.data[:], 
                                         g_ref, 
                                         decimal=decimal)

if __name__ == "__main__":
    test_multiorbital_F_wk()
