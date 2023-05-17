# ----------------------------------------------------------------------

import numpy as np

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import eliashberg_g_delta_g_product
from triqs_tprf.lattice import dlr_on_imfreq

from triqs.gf import Gf, MeshImFreq, MeshBrillouinZone
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRCoeffs
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs.gf.gf_fnt import dlr_coeffs_from_dlr_imfreq

# ----------------------------------------------------------------------

def compare_g_Dwk_and_g_wk(g_Dwk, g_wk, decimal=7):
    wmesh, kmesh = g_wk.mesh.components
    DLRwmesh = g_Dwk.mesh.components[0]

    g_ref_wk =Gf(mesh=g_wk.mesh, target_shape=g_wk.target_shape)
    for k in kmesh:
        g_Dw = Gf(mesh=DLRwmesh, target_shape=g_Dwk.target_shape)
        g_Dw.data[:] = g_Dwk.data[:,k.linear_index,:]
        g_Dc = dlr_coeffs_from_dlr_imfreq(g_Dw)
        g_ref_wk[:,k] = dlr_on_imfreq(g_Dc, wmesh)
    
    np.testing.assert_array_almost_equal(g_wk.data[:], g_ref_wk.data[:], decimal=decimal)

def test_g_delta_g_product_dlr():
    """ Some test description
    Author: Yann in 't Veld (2023) """ 
    
    mu = 0.5
    g2 = 0.4
    wD = 0.1
    beta = 10.0
    nk = 10

    nw = 50
    lamb = 100
    eps = 1e-10

    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrillouinZone(bz, np.diag(np.array([nk, nk, nk], dtype=int)))
    
    wmesh = MeshImFreq(beta, 'Fermion', nw)
    DLRwmesh = MeshDLRImFreq(beta, 'Fermion', lamb, eps)

    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        Enk.data[k.linear_index,:] = 0.1 * knorm**2.0

    print(np.max(Enk.data[:].real))
    print(np.max(Enk.data[:].real) * beta)

    g0_wk = lattice_dyson_g0_wk(mu, Enk, wmesh)
    g0_Dwk = lattice_dyson_g0_wk(mu, Enk, DLRwmesh)
    compare_g_Dwk_and_g_wk(g0_Dwk, g0_wk)

    print('--> delta_wk')
    delta_wk = Gf(mesh=g0_wk.mesh, target_shape=g0_wk.target_shape)
    for w in wmesh:
        wii = w.linear_index
        delta_wk.data[wii,:] = 1.0 / (w.value + Enk.data[:])

    delta_Dwk = Gf(mesh=g0_Dwk.mesh, target_shape=g0_Dwk.target_shape)
    for w in DLRwmesh:
        wii = w.linear_index
        delta_Dwk.data[wii,:] = 1.0 / (w.value + Enk.data[:])

    compare_g_Dwk_and_g_wk(delta_Dwk, delta_wk)

    print('--> eliashberg_g_delta_g_product')
    F_wk = eliashberg_g_delta_g_product(g0_wk, delta_wk)
    F_Dwk = eliashberg_g_delta_g_product(g0_Dwk, delta_Dwk)

    print(F_wk.data[nw,0:20,0,0])

    compare_g_Dwk_and_g_wk(F_Dwk, F_wk)

if __name__ == "__main__":
    test_g_delta_g_product_dlr()
