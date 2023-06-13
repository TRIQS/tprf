# ----------------------------------------------------------------------

import numpy as np

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import eliashberg_g_delta_g_product, dynamic_and_constant_to_tr
from triqs_tprf.lattice import eliashberg_product_fft, eliashberg_product_fft_constant
from triqs_tprf.lattice import dlr_on_imfreq, dlr_on_imtime

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr


from triqs.gf import Gf, MeshImFreq, MeshBrZone
from triqs.gf.meshes import MeshDLRImFreq, MeshDLR
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs.gf.gf_factories import make_gf_dlr

# ----------------------------------------------------------------------

def compare_g_Dwk_and_g_wk(g_Dwk, g_wk, decimal=7):
    wmesh, kmesh = g_wk.mesh.components
    DLRwmesh = g_Dwk.mesh.components[0]

    g_ref_wk =Gf(mesh=g_wk.mesh, target_shape=g_wk.target_shape)
    for k in kmesh:
        g_Dw = Gf(mesh=DLRwmesh, target_shape=g_Dwk.target_shape)
        g_Dw.data[:] = g_Dwk.data[:,k.data_index,:]
        g_Dc = make_gf_dlr(g_Dw)
        g_ref_wk[:,k] = dlr_on_imfreq(g_Dc, wmesh)
    
    np.testing.assert_array_almost_equal(g_wk.data[:], g_ref_wk.data[:], decimal=decimal)

def compare_g_Dtk_and_g_tk(g_Dtk, g_tk, decimal=7):
    tmesh, kmesh = g_tk.mesh.components
    DLRtmesh = g_Dtk.mesh.components[0]

    print(kmesh)
    for k in kmesh:
        print(k)
    
    g_ref_tk =Gf(mesh=g_tk.mesh, target_shape=g_tk.target_shape)
    for k in kmesh:
        g_Dt = Gf(mesh=DLRtmesh, target_shape=g_Dtk.target_shape)
        g_Dt.data[:] = g_Dtk.data[:,k.data_index,:]
        g_Dc = make_gf_dlr(g_Dt)
        g_ref_tk[:,k] = dlr_on_imtime(g_Dc, tmesh)
    
    np.testing.assert_array_almost_equal(g_tk.data[:], g_ref_tk.data[:], decimal=decimal)


def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)


def eliashberg_compare_dlr_and_direct():
    """ Some test description
    Author: Yann in 't Veld (2023) """ 
    
    mu = 0.5
    g2 = 0.4
    wD = 0.1
    beta = 2.0
    nk = 3

    nw = 50
    lamb = 10.
    eps = 1e-8

    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, nk])
    
    wmesh = MeshImFreq(beta, 'Fermion', nw)
    DLRwmesh = MeshDLRImFreq(beta, 'Fermion', lamb, eps)

    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        Enk.data[k.data_index,:] = 0.1 * knorm**2.0

    print(np.max(Enk.data[:].real))
    print(np.max(Enk.data[:].real) * beta)

    g0_wk = lattice_dyson_g0_wk(mu, Enk, wmesh)
    g0_Dwk = lattice_dyson_g0_wk(mu, Enk, DLRwmesh)
    compare_g_Dwk_and_g_wk(g0_Dwk, g0_wk)

    print('--> delta_wk')
    delta_wk = Gf(mesh=g0_wk.mesh, target_shape=g0_wk.target_shape)
    for w in wmesh:
        wii = w.data_index
        delta_wk.data[wii,:] = 1.0 / (w.value + Enk.data[:])

    delta_Dwk = Gf(mesh=g0_Dwk.mesh, target_shape=g0_Dwk.target_shape)
    for w in DLRwmesh:
        wii = w.data_index
        delta_Dwk.data[wii,:] = 1.0 / (w.value + Enk.data[:])

    compare_g_Dwk_and_g_wk(delta_Dwk, delta_wk)

    print('--> eliashberg_g_delta_g_product')
    F_wk = eliashberg_g_delta_g_product(g0_wk, delta_wk)
    F_Dwk = eliashberg_g_delta_g_product(g0_Dwk, delta_Dwk)

    print(F_wk.data[nw,0:20,0,0])

    compare_g_Dwk_and_g_wk(F_Dwk, F_wk)

    print('--> setup interaction vertex')
    numesh = MeshImFreq(beta, 'Boson', nw)
    DLRnumesh = MeshDLRImFreq(beta, 'Boson', lamb, eps)

    I_wk = Gf(mesh=MeshProduct(numesh, kmesh), target_shape=[1]*4)
    for nu in numesh:
        nuii = nu.data_index
        I_wk.data[nuii,:] = ElectronPhononInteraction(nu.value, g2 ,wD)

    I_Dwk = Gf(mesh=MeshProduct(DLRnumesh, kmesh), target_shape=[1]*4)
    for nu in DLRnumesh:
        nuii = nu.data_index
        I_Dwk.data[nuii,:] = ElectronPhononInteraction(nu.value, g2 ,wD)

    I_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_k.data[:] = 0.0
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        if(np.isclose(knorm, 0.0)): break
        I_k.data[:] = 1.0 / knorm 

    compare_g_Dwk_and_g_wk(I_Dwk, I_wk)

    print("--> dynamic_and_constant_to_tr")
    I_dyn_tr, I_r_ref = dynamic_and_constant_to_tr(I_wk, I_k)
    I_dyn_Dtr, I_r = dynamic_and_constant_to_tr(I_Dwk, I_k)

    compare_g_Dtk_and_g_tk(I_dyn_Dtr, I_dyn_tr, decimal=6)
    np.testing.assert_array_almost_equal(I_r.data[:], I_r_ref.data[:])

    print("--> eliashberg_product_fft")
    F_wr = fourier_wk_to_wr(F_wk)
    F_tr = fourier_wr_to_tr(F_wr)

    F_Dwr = fourier_wk_to_wr(F_Dwk)
    F_Dtr = fourier_wr_to_tr(F_Dwr)

    delta_wk_out = eliashberg_product_fft(I_dyn_tr, I_r, g0_wk, delta_wk)
    delta_Dwk_out = eliashberg_product_fft(I_dyn_Dtr, I_r, g0_Dwk, delta_Dwk)
    compare_g_Dwk_and_g_wk(delta_Dwk_out, delta_wk_out)


    print("--> eliashberg_product_fft_constant")

    delta_wk_out_const = eliashberg_product_fft_constant(I_r, g0_wk, delta_wk)
    delta_Dwk_out_const = eliashberg_product_fft_constant(I_r, g0_Dwk, delta_Dwk)
    compare_g_Dwk_and_g_wk(delta_Dwk_out_const, delta_wk_out_const)

if __name__ == "__main__":
    eliashberg_compare_dlr_and_direct()
