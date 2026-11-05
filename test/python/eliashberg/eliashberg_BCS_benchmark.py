import numpy as np
import sys
import time
import os

import matplotlib.pyplot as plt

from triqs.gfs import *
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk

from triqs_tprf.eliashberg import solve_eliashberg, semi_random_initial_delta
from triqs_tprf.eliashberg import solve_eliashberg_nonlinear

def BCS_gap_equation(E, U, beta, delta=0.0):
    # Additional factor of 0.5 to account for the same factor
    # in the Eliashberg functionality of TPRF
    E_gapped = np.sqrt(E**2.0 + delta**2.0)
    out = -0.5 * U * np.tanh( 0.5 * beta * E_gapped ) / (2.0 * E_gapped)
    return out

def test_compare_Eliashberg_to_BCS():
    nk = 1  # Use a completely local model, so nk=1 is enough
    beta = 100.0
    dlr_wmax = 50.0
    dlr_eps = 1e-12
    E = -1.0
    U = -5.0

    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, nk, nk])
    
    wmesh_f  = MeshDLRImFreq(beta, 'Fermion', dlr_wmax, dlr_eps)
    wmesh_b = MeshDLRImFreq(beta, 'Boson', dlr_wmax, dlr_eps)

    print('--> construct local model')
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    Enk.data[:] = E

    g0_wk = lattice_dyson_g0_wk(mu=0.0, e_k=Enk, mesh=wmesh_f)

    I_wk = Gf(mesh=MeshProduct(wmesh_b, kmesh), target_shape=[1]*4)
    I_wk.data[:] = U
    
    I_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_k.data[:] = U

    print('--> linearized Eliashberg')
    delta0_wk = g0_wk.copy()
    delta0_wk.data[:] = 0.80

    vals, vecs = solve_eliashberg(I_wk, 
                                  g0_wk, 
                                  initial_delta=delta0_wk, 
                                  Gamma_pp_const_k=I_k, 
                                  product="FFT", 
                                  solver="IRAM", 
                                  tol=1e-4)
    leadingIndex = np.argmax(np.real(vals))
    leadingEv = vals[leadingIndex]
    delta_lin_wk = vecs[leadingIndex]

    assert len(np.unique(delta_lin_wk.data[:])) == 1, "Linearized gap was not constant!"

    leadingEv_BCS = BCS_gap_equation(E, U, beta, delta=0.0)

    print("  Leading eigenvalues:")
    print("  Eliashberg:", leadingEv)
    print("  BCS:", leadingEv_BCS)

    assert np.isclose(leadingEv, leadingEv_BCS, 1e-2), "Eliashberg and BCS did not agree!"

    print('--> non-linearized Eliashberg')
    delta_nonlin_wk = solve_eliashberg_nonlinear(I_wk, 
                                                 g0_wk, 
                                                 initial_delta=delta0_wk, 
                                                 Gamma_pp_const_k=I_k, 
                                                 product="FFT", 
                                                 steffensens_method=True, 
                                                 tol=1e-4,
                                                 max_it=100)

    assert len(np.unique(delta_nonlin_wk.data[:])) == 1, "Gap was not constant!"
    delta = np.unique(delta_nonlin_wk.data[:,0,0])[0]

    print("  delta =", delta)

    LHS_BCS = BCS_gap_equation(E, U, beta, delta)
    print("  LHS BCS gap equation =", LHS_BCS)

    assert np.isclose(LHS_BCS, 1.0, 1e-2), "Full Eliashberg and BCS did not agree!"



if __name__ == "__main__":
    test_compare_Eliashberg_to_BCS()
