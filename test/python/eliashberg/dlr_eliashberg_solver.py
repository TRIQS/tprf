# ----------------------------------------------------------------------

import numpy as np

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import dlr_on_imfreq, dlr_on_imtime
from triqs_tprf.eliashberg import solve_eliashberg, semi_random_initial_delta

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr


from triqs.gf import Gf, MeshImFreq, MeshBrZone
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRCoeffs
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs.gf.gf_fnt import dlr_coeffs_from_dlr_imfreq, dlr_coeffs_from_dlr_imtime

# ----------------------------------------------------------------------

def dlr_wk_on_imfreq_wk(g_Dwk, wmesh):
    DLRwmesh, kmesh = g_Dwk.mesh.components
    g_wk =Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=g_Dwk.target_shape)
    g_wk.data[:] = 0.0

    for k in kmesh:
        g_Dw = Gf(mesh=DLRwmesh, target_shape=g_Dwk.target_shape)
        g_Dw.data[:] = g_Dwk.data[:,k.data_index,:]
        g_Dc = dlr_coeffs_from_dlr_imfreq(g_Dw)
        g_wk[:,k] = dlr_on_imfreq(g_Dc, wmesh)

    return g_wk


def compare_g_Dwk_and_g_wk(g_Dwk, g_wk, decimal=7):
    wmesh, kmesh = g_wk.mesh.components
    DLRwmesh = g_Dwk.mesh.components[0]

    g_ref_wk = dlr_wk_on_imfreq_wk(g_Dwk, wmesh)

    np.testing.assert_array_almost_equal(g_wk.data[:], g_ref_wk.data[:], decimal=decimal)

def compare_g_Dtk_and_g_tk(g_Dtk, g_tk, decimal=7):
    tmesh, kmesh = g_tk.mesh.components
    DLRtmesh = g_Dtk.mesh.components[0]

    g_ref_tk =Gf(mesh=g_tk.mesh, target_shape=g_tk.target_shape)
    for k in kmesh:
        g_Dt = Gf(mesh=DLRtmesh, target_shape=g_Dtk.target_shape)
        g_Dt.data[:] = g_Dtk.data[:,k.data_index,:]
        g_Dc = dlr_coeffs_from_dlr_imtime(g_Dt)
        g_ref_tk[:,k] = dlr_on_imtime(g_Dc, tmesh)
    
    np.testing.assert_array_almost_equal(g_tk.data[:], g_ref_tk.data[:], decimal=decimal)


def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)


def test_dlr_eliashberg_solver():
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
    kmesh = MeshBrZone(bz, np.diag(np.array([nk, nk, nk], dtype=int)))
    
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
    delta0_Dwk = semi_random_initial_delta(g0_Dwk)
    delta0_wk = dlr_wk_on_imfreq_wk(delta0_Dwk, wmesh)

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

    print("--> solve_eliashberg using DLR")
    vals_dlr, vecs_dlr = solve_eliashberg(I_Dwk, g0_Dwk, initial_delta=delta0_Dwk, Gamma_pp_const_k=I_k, product="FFT", solver="IRAM")
    leadingIndex_dlr = np.argmax(np.real(vals_dlr))
    delta_Dwk_out = vecs_dlr[leadingIndex_dlr]
    delta_wk_out_dlr = dlr_wk_on_imfreq_wk(delta_Dwk_out, wmesh)
    delta_wk_out_dlr.data[:] /= delta_wk_out_dlr.data[len(wmesh)//2,0,0,0]

    print("--> solve_eliashberg directly")
    vals, vecs = solve_eliashberg(I_wk, g0_wk, initial_delta=delta0_wk, Gamma_pp_const_k=I_k, product="FFT", solver="IRAM")
    leadingIndex = np.argmax(np.real(vals))
    delta_wk_out = vecs[leadingIndex]
    delta_wk_out.data[:] /= delta_wk_out.data[len(wmesh)//2,0,0,0]
    print("")

    print("--> compare DLR and direct")
    print("Leading eigenvalue DLR:   ", vals_dlr[leadingIndex_dlr])
    print("Leading eigenvalue direct:", vals[leadingIndex])
    np.testing.assert_array_almost_equal(np.sort(vals_dlr), np.sort(vals))
    np.testing.assert_array_almost_equal(delta_wk_out_dlr.data[:], delta_wk_out.data[:])

if __name__ == "__main__":
    test_dlr_eliashberg_solver()
