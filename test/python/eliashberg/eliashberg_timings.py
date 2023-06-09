# ----------------------------------------------------------------------

import numpy as np
import time

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import eliashberg_g_delta_g_product, dynamic_and_constant_to_tr
from triqs_tprf.lattice import eliashberg_product_fft, eliashberg_product_fft_constant

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr

from triqs.gf import Gf, MeshImFreq, MeshBrillouinZone
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRCoeffs
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)


def eliashberg_timings():
    """ Some test description
    Author: Yann in 't Veld (2023) """ 
    
    mu = 0.5
    g2 = 0.4
    wD = 0.1
    beta = 10.0
    nk = 10

    lamb = 1000
    eps = 1e-10
    nw = int(lamb / (2.0 * np.pi) - 0.5) #50

    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrillouinZone(bz, np.diag(np.array([nk, nk, nk], dtype=int)))
    
    wmesh = MeshImFreq(beta, 'Fermion', nw)
    DLRwmesh = MeshDLRImFreq(beta, 'Fermion', lamb, eps)

    print("  size linear mesh: %i"%len(wmesh))
    print("  size DLR mesh:    %i"%len(DLRwmesh))



    print('--> lattice_dyson_g0_wk')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        Enk.data[k.linear_index,:] = 0.1 * knorm**2.0

    g0_wk = lattice_dyson_g0_wk(mu, Enk, wmesh)
    g0_Dwk = lattice_dyson_g0_wk(mu, Enk, DLRwmesh)



    print('--> delta_wk')
    delta_wk = Gf(mesh=g0_wk.mesh, target_shape=g0_wk.target_shape)
    for w in wmesh:
        wii = w.linear_index
        delta_wk.data[wii,:] = 1.0 / (w.value + Enk.data[:])

    delta_Dwk = Gf(mesh=g0_Dwk.mesh, target_shape=g0_Dwk.target_shape)
    for w in DLRwmesh:
        wii = w.linear_index
        delta_Dwk.data[wii,:] = 1.0 / (w.value + Enk.data[:])




    print('--> eliashberg_g_delta_g_product')
    start = time.time()
    F_wk = eliashberg_g_delta_g_product(g0_wk, delta_wk)
    end = time.time()
    print("  Linear: %6.4f secs"%(end-start))

    start = time.time()
    F_Dwk = eliashberg_g_delta_g_product(g0_Dwk, delta_Dwk)
    end = time.time()
    print("  DLR:    %6.4f secs"%(end-start))





    print('--> setup interaction vertex')
    numesh = MeshImFreq(beta, 'Boson', nw)
    DLRnumesh = MeshDLRImFreq(beta, 'Boson', lamb, eps)

    I_wk = Gf(mesh=MeshProduct(numesh, kmesh), target_shape=[1]*4)
    for nu in numesh:
        nuii = nu.linear_index
        I_wk.data[nuii,:] = ElectronPhononInteraction(nu.value, g2 ,wD)

    I_Dwk = Gf(mesh=MeshProduct(DLRnumesh, kmesh), target_shape=[1]*4)
    for nu in DLRnumesh:
        nuii = nu.linear_index
        I_Dwk.data[nuii,:] = ElectronPhononInteraction(nu.value, g2 ,wD)

    I_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_k.data[:] = 0.0
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        if(np.isclose(knorm, 0.0)): break
        I_k.data[:] = 1.0 / knorm 


    print("--> dynamic_and_constant_to_tr")
    I_dyn_tr, I_r_ref = dynamic_and_constant_to_tr(I_wk, I_k)
    I_dyn_Dtr, I_r = dynamic_and_constant_to_tr(I_Dwk, I_k)




    print("--> eliashberg_product_fft")
    F_wr = fourier_wk_to_wr(F_wk)
    F_tr = fourier_wr_to_tr(F_wr)

    F_Dwr = fourier_wk_to_wr(F_Dwk)
    F_Dtr = fourier_wr_to_tr(F_Dwr)

    start = time.time()
    delta_wk_out = eliashberg_product_fft(I_dyn_tr, I_r, g0_wk, delta_wk)
    end = time.time()
    print("  Linear: %6.4f secs"%(end-start))

    start = time.time()
    delta_Dwk_out = eliashberg_product_fft(I_dyn_Dtr, I_r, g0_Dwk, delta_Dwk)
    end = time.time()
    print("  DLR:    %6.4f secs"%(end-start))
    


    print("--> eliashberg_product_fft_constant")

    start = time.time()
    delta_wk_out_const = eliashberg_product_fft_constant(I_r, g0_wk, delta_wk)
    end = time.time()
    print("  Linear: %6.4f secs"%(end-start))

    start = time.time()
    delta_Dwk_out_const = eliashberg_product_fft_constant(I_r, g0_Dwk, delta_Dwk)
    end = time.time()
    print("  DLR:    %6.4f secs"%(end-start))


if __name__ == "__main__":
    eliashberg_timings()
