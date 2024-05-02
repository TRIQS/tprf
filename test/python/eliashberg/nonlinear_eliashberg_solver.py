import numpy as np
import sys
import time
import os

import matplotlib.pyplot as plt

from triqs.gf import *
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk, dynamic_and_constant_to_tr

from triqs_tprf.eliashberg import solve_eliashberg, semi_random_initial_delta
from triqs_tprf.eliashberg import solve_eliashberg_nonlinear

kB = 8.6173427909 * 10**(-5.0)

def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)

def SolveEliashbergAtTemperature(T):
    print("T = %.5f K"%T)
    
    mu = 1.0
    g2 = 0.3
    wD = 0.3
    nk = 8

    wmax = 20.
    eps = 1e-10
    
    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, nk])
    
    beta = 1.0 / (kB * T)
    wmesh_f  = MeshDLRImFreq(beta, 'Fermion', wmax, eps)
    wmesh_b = MeshDLRImFreq(beta, 'Boson', wmax, eps)
    
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        Enk.data[k.data_index,:] = 0.1 * knorm**2.0

    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=Enk, mesh=wmesh_f)

    print('--> setup interaction vertex')
    I_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_k.data[:] = 0.0
    for k in kmesh:
        knorm = np.linalg.norm(k.value)
        if(np.isclose(knorm, 0.0)): continue
        I_k.data[:] = 0.0 #1.0 / knorm

    I_wk = Gf(mesh=MeshProduct(wmesh_b, kmesh), target_shape=[1]*4)
    for w in wmesh_b:
        wii = w.data_index
        I_wk.data[wii,:] = ElectronPhononInteraction(w.value, g2 ,wD) + I_k.data[:]
        
    print('--> delta_wk')
    delta0_wk = semi_random_initial_delta(g0_wk)
        
    print("--> solve_eliashberg")
    vals, vecs = solve_eliashberg(I_wk, g0_wk, initial_delta=delta0_wk, Gamma_pp_const_k=I_k, product="FFT", solver="IRAM")
    leadingIndex = np.argmax(np.real(vals))
    leadingEv = vals[leadingIndex]
    delta_lin_wk = vecs[leadingIndex]

    print("--> solve_eliashberg_nonlinear")
    delta_nonlin_wk = solve_eliashberg_nonlinear(I_wk, g0_wk, initial_delta=delta0_wk, Gamma_pp_const_k=I_k, product="FFT")

    return leadingEv, delta_lin_wk, delta_nonlin_wk


def GetDeltaAtiw0(delta_wk, T):
    k0ind = Idx([0,0,0])
    delta_c = make_gf_dlr(delta_wk[:,k0ind])
    
    beta = 1.0 / (kB * T)
    return delta_c(MatsubaraFreq(0, beta, "Fermion"))[0,0].real


def test_nonlinear_Eliashberg_solver(verbose=False):
    """ Test the non-linear Eliahsberg solver against the linearized Eliashberg solver
    Author: Yann in 't Veld (2024) """

    if verbose:
        TLst = [200, 225, 250, 260, 300]
    else:
        TLst = [200, 300]
    datLst = np.array([SolveEliashbergAtTemperature(T) for T in TLst])
    
    leadingEvLst = datLst[:,0]
    delta_lin_wkLst = datLst[:,1]
    delta_nonlin_wkLst = datLst[:,2]
    
    delta_iw0Lst = [GetDeltaAtiw0(delta_wk, TLst[ii]) for ii,delta_wk in enumerate(delta_nonlin_wkLst)]
    
    for ii in range(len(leadingEvLst)):
        ev = leadingEvLst[ii]
        delta_iw0 = delta_iw0Lst[ii]
        
        if ev > 1.0:
            assert delta_iw0 > 1e-6
        if ev < 1.0:
            assert np.isclose(delta_iw0, 0.0, atol=1e-6)
    
    if verbose:
        leadingEvLst = np.array(leadingEvLst, dtype=float)
        Tc = np.interp(1.0, leadingEvLst[::-1], TLst[::-1])
        print("Tc =", Tc, "K")
        
        fig = plt.figure(figsize=(5.0,3.0),dpi=300)

        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()

        ax1.plot(TLst, leadingEvLst, marker='.', color='cornflowerblue')
        ax1.axhline(y=1.0, linewidth=0.5, color='gray', linestyle=':')
        ax1.set_ylabel("$\lambda$", color="cornflowerblue")

        ax2.plot(TLst, delta_iw0Lst, marker='.', color='firebrick')
        ax2.axhline(y=0.0, linewidth=0.5, color='black', linestyle='--')
        ax2.set_ylabel("$\Delta_{k=0}(i\omega_n=0)$", color="firebrick")

        ax1.set_xlabel("T (K)")

        plt.tight_layout()
        plt.show()
    
    
    
if __name__ == "__main__":
    test_nonlinear_Eliashberg_solver(verbose=False)
