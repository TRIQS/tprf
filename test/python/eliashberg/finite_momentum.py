import numpy as np
import sys
import time
import os

import matplotlib.pyplot as plt

from triqs.gfs import *
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk, dynamic_and_constant_to_tr

from triqs_tprf.eliashberg import solve_eliashberg, semi_random_initial_delta
from triqs_tprf.eliashberg import solve_eliashberg_nonlinear

kB = 8.6173427909 * 10**(-5.0)

def AnalyticalGapEquation(Enk, Hnkmq, mu, beta, U, delta_k=None):
    """The gap equation for a constant interaction U.
    Here the Matsubara sum was done analytically.
    The momentum sum is done numerically.
    For a correct solution delta_k this function is expected to return unity."""

    if delta_k is None:
        delta_k = Enk.copy()
        delta_k.data[:] = 0.0

    def poles(Enk, Hnk, mu, delta_k, sign=1.0):
        E = Enk.data[:,0,0] - mu
        H = Hnk.data[:,0,0] + mu
        D = delta_k.data[:,0,0]

        sqrtval = (E + H)**2.0\
                + 4.0 * (-E*H + np.abs(D)**2.0)
        return 0.5 * (E + H + sign * np.sqrt(sqrtval))

    kmesh = Enk.mesh
    fermi = lambda E: 0.5 * (1.0 - np.tanh(0.5 * beta * E))
    
    zpos = poles(Enk, Hnkmq, mu, delta_k, sign=1.0)
    zneg = poles(Enk, Hnkmq, mu, delta_k, sign=-1.0)
    integrant = (1.0 - fermi(zpos) - fermi(-zneg)) / (zpos - zneg)
    
    integval = np.sum(integrant) / len(kmesh)
    funcval = - 0.5 * U * integval
    return funcval



def backfold_k(kvec, kmesh):
    """Backfold a k-vector into the first Brillouin Zone

    Parameters
    ----------
    kvec : numpy.ndarray [shape=(3)]
       The k-vector to backfold into the first Brillouin Zone
    kmesh: MeshBrZone
       The Brillouin Zone mesh on which to backfold

    Returns
    -------
    kvecsFolded: numpy.ndarray [shape=(3)]
        The folded k-vector within the first Brillouin Zone
    """

    # get the k-vector in internal units
    kvecInt = kvec @ np.linalg.inv(kmesh.bz.units)

    # backfold the k-vector in internal units
    # to [-0.5, 0.5] for each dimension
    for ii in range(len(kvecInt)):
        kvecInt[ii] = (kvecInt[ii] + 0.5)%1.0 - 0.5

    # get backfolded k-vector in physical units 1/ang
    kvecFolded = kvecInt @ kmesh.bz.units
    return kvecFolded

def SolveEliashbergAtTemperature(T, q_fmp):
    print("T = %.5f K"%T)
    print("Q = %.5f"%q_fmp)
    
    mu = 1.0
    nk = 8
    U = -1.0
    wmax = 50.
    eps = 1e-10
    
    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, 1, 1])
    
    beta = 1.0 / (kB * T)
    wmesh_f  = MeshDLRImFreq(beta, 'Fermion', wmax, eps)
    wmesh_b = MeshDLRImFreq(beta, 'Boson', wmax, eps)
    
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)   # Electron Gf
    Hnkmq = Gf(mesh=kmesh, target_shape=[1]*2) # Hole Gf
    for k in kmesh:
        kvec = backfold_k(k.value, kmesh)
        knorm = np.linalg.norm(kvec)
        Enk.data[k.data_index,:] = 0.1 * knorm**2.0

        kvec = backfold_k(k.value - q_fmp, kmesh)
        knorm = np.linalg.norm(kvec)
        Hnkmq.data[k.data_index,:] = -0.1 * knorm**2.0

    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=Enk, mesh=wmesh_f)
    g0bar_wk = lattice_dyson_g0_wk(mu=-mu, e_k=Hnkmq, mesh=wmesh_f)

    print('--> setup interaction vertex')
    I_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_k.data[:] = 0.0
    for k in kmesh:
        I_k.data[:] = U

    I_wk = Gf(mesh=MeshProduct(wmesh_b, kmesh), target_shape=[1]*4)
    I_wk.data[:] = U
        
    print('--> delta_wk')
    delta0_wk = semi_random_initial_delta(g0_wk)
        
    print("--> solve_eliashberg")
    vals, vecs = solve_eliashberg(I_wk, 
                                  g0_wk, 
                                  initial_delta=delta0_wk, 
                                  Gamma_pp_const_k=I_k, 
                                  product="FFT", 
                                  solver="IRAM", 
                                  tol=1e-8,
                                  gbar_wk=g0bar_wk)
    leadingIndex = np.argmax(np.real(vals))
    leadingEv = vals[leadingIndex]
    delta_lin_wk = vecs[leadingIndex]
    print("  eigenvalue:", leadingEv)

    leadingEv_ana = AnalyticalGapEquation(Enk, Hnkmq, mu, beta, U)
    print("  analytical eigenvalue:", leadingEv_ana)

    assert np.isclose(leadingEv, leadingEv_ana, atol=1e-2), "Eigenvalue is not as expected"

    print("--> solve_eliashberg_nonlinear")
    delta_nonlin_wk = solve_eliashberg_nonlinear(I_wk, 
                                                 g0_wk, 
                                                 initial_delta=delta_lin_wk, 
                                                 Gamma_pp_const_k=I_k, 
                                                 product="FFT", 
                                                 tol=1e-8,
                                                 gbar_wk=g0bar_wk)
   
    delta_iw0 = np.abs( GetDeltaAtiw0(delta_nonlin_wk, T) )
    print("  delta(iw=0):", delta_iw0)

    if leadingEv > 1.0:
        assert delta_iw0 > 1e-6, "Gap should be finite below Tc!"

        delta_k = Enk.copy()
        delta_k.data[:] = delta_nonlin_wk.data[0,:]
        LHS_ana = AnalyticalGapEquation(Enk, Hnkmq, mu, beta, U, delta_k)
        print("  Analytical gap LHS:", LHS_ana)

        assert np.isclose(LHS_ana, 1.0, atol=1e-2), "Gap is not as expected"

    if leadingEv < 1.0:
        assert np.isclose(delta_iw0, 0.0, atol=1e-6), "Gap should vanish above Tc!"

    return leadingEv, delta_lin_wk, delta_nonlin_wk


def GetDeltaAtiw0(delta_wk, T):
    k0ind = Idx([0,0,0])
    delta_c = make_gf_dlr(delta_wk[:,k0ind])
    
    beta = 1.0 / (kB * T)
    return delta_c(MatsubaraFreq(0, beta, "Fermion"))[0,0].real


def test_nonlinear_Eliashberg_solver(verbose=False):
    """ Test the linearized and non-linearized Eliashberg solvers against
    an expression where the Matsubara sums are done analytically.
    This is done for a local interaction U.
    We also include a finite momentum pairing (FMP) to the Cooper pairs,
    by shifting the time-reversed Green's function.
    See also: https://doi.org/10.1038/s41535-024-00706-7
    Author: Yann in 't Veld (2026) """

    T = 100.0
    if verbose:
        QLst = [0.00, 0.03, 0.06, 0.09, 0.12, 0.15]
    else:
        QLst = [0.06, 0.12]
    datLst = np.array([SolveEliashbergAtTemperature(T,Q) for Q in QLst])
    
    leadingEvLst = datLst[:,0]
    delta_lin_wkLst = datLst[:,1]
    delta_nonlin_wkLst = datLst[:,2]
    
    delta_iw0Lst = [GetDeltaAtiw0(delta_wk, T) for ii,delta_wk in enumerate(delta_nonlin_wkLst)]
    
    for ii in range(len(leadingEvLst)):
        ev = leadingEvLst[ii]
        delta_iw0 = np.abs(delta_iw0Lst[ii])
        
    
    if verbose:
        leadingEvLst = np.array(leadingEvLst, dtype=float)
        Qc = np.interp(1.0, leadingEvLst[::-1], QLst[::-1])
        print("Qc =", Qc)
        
        fig = plt.figure(figsize=(5.0,3.0),dpi=300)

        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()

        ax1.plot(QLst, leadingEvLst, marker='.', color='cornflowerblue')
        ax1.axhline(y=1.0, linewidth=0.5, color='gray', linestyle=':')
        ax1.set_ylabel("$\lambda$", color="cornflowerblue")

        ax2.plot(QLst, np.abs(delta_iw0Lst), marker='.', color='firebrick')
        ax2.axhline(y=0.0, linewidth=0.5, color='black', linestyle='--')
        ax2.set_ylabel("$|\Delta_{k=0}(i\omega_n=0)$|", color="firebrick")

        ax1.set_xlabel("Q")

        plt.tight_layout()
        plt.show()
    
    
    
if __name__ == "__main__":
    test_nonlinear_Eliashberg_solver(verbose=False)
