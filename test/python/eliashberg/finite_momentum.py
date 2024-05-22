import numpy as np

from triqs.gf import *
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.eliashberg import solve_eliashberg, semi_random_initial_delta
from triqs_tprf.eliashberg import solve_eliashberg_nonlinear

kB = 8.6173427909 * 10**(-5.0)


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

def Xi(kvec, mu):
    kmax = np.pi
    knorm = kvec.copy()
    knorm[kvec >= kmax] = knorm[kvec >= kmax] - 2.0*kmax
    knorm[kvec < -kmax] = knorm[kvec < -kmax] + 2.0*kmax
    return knorm**2.0 - mu

def poles(kvec, phi, qvec, mu, sign=1.0):
    xi_kpq = Xi(kvec+0.5*qvec, mu)
    xi_mkpq = Xi(-kvec+0.5*qvec, mu)
    
    sqrtval = (xi_kpq - xi_mkpq)**2.0 + 4.0 * (xi_kpq * xi_mkpq + np.abs(phi)**2.0)
    return 0.5 * (xi_kpq - xi_mkpq + sign * np.sqrt(sqrtval))

def AnalyticalMatsubaraSumGapEquation(phi, beta, qvec, U, mu):
    """The gap equation for a constant interaction U.
    Here the Matsubara sum was done analytically. The momentum sum is done using 
    the Trapezoid method implemented in Numpy.
    For a correct solution phi this function is expected to return unity."""
    
    kmax = np.pi
    kmesh = np.linspace(-kmax, kmax, 500)
    fermi = lambda E: 0.5 * (1.0 - np.tanh(0.5 * beta * E))
    
    zpos = poles(kmesh, phi, qvec, mu, sign=1.0)
    zneg = poles(kmesh, phi, qvec, mu, sign=-1.0)
    integrant = (1.0 - fermi(zpos) - fermi(-zneg)) / (zpos - zneg)
    
    integval = np.trapz(integrant, kmesh) / (2.0 * np.pi)
    funcval = - 0.5 * U * integval
    return funcval

def SolveEliashbergAtTemperature(T, qind):
    print("T = %.5f K"%T)

    mu = 2.0
    U = -3.5
    nk = 500

    wmax = 100.
    eps = 1e-8

    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, 1, 1])
    q_fmp = np.linalg.norm(kmesh[qind].value)

    beta = 1.0 / (kB * T)
    wmesh_f  = MeshDLRImFreq(beta, 'Fermion', wmax, eps)
    wmesh_b = MeshDLRImFreq(beta, 'Boson', wmax, eps)

    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        kfolded = backfold_k(k.value, kmesh)
        knorm = np.linalg.norm(kfolded)
        Enk.data[k.data_index,:] = knorm**2.0

    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=Enk, mesh=wmesh_f)

    print('--> setup interaction vertex')
    I_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_k.data[:] = U

    I_wk = Gf(mesh=MeshProduct(wmesh_b, kmesh), target_shape=[1]*4)
    I_wk.data[:] = U

    print('--> delta_wk')
    delta0_wk = semi_random_initial_delta(g0_wk)

    print("--> solve_eliashberg")
    vals, vecs = solve_eliashberg(I_wk, g0_wk, initial_delta=delta0_wk, Gamma_pp_const_k=I_k, product="FFT", solver="IRAM", tol=1e-6, fmp_index=qind)
    leadingIndex = np.argmax(np.real(vals))
    leadingEv = vals[leadingIndex]
    delta_lin_wk = vecs[leadingIndex]
    
    # Expect constant eigenvector
    delta_lin = np.unique( np.round( delta_lin_wk.data[:], decimals=6) )
    assert np.shape(delta_lin) == (1,)
    delta_lin = np.abs(delta_lin[0])
    
    # Compare with leading eigenvalue evaluated from analytical Matsubara summations
    leadingEv_ana = AnalyticalMatsubaraSumGapEquation(0.0, beta, q_fmp, U, mu)
    np.testing.assert_array_almost_equal(leadingEv, leadingEv_ana, decimal=3)

    print("--> solve_eliashberg_nonlinear")
    delta_nonlin_wk = solve_eliashberg_nonlinear(I_wk, g0_wk, initial_delta=delta0_wk, Gamma_pp_const_k=I_k, product="FFT", tol=1e-6, mixing_frac=0.2, fmp_index=qind)
                
    # Expect constant eigenvector
    delta_nonlin = np.unique( np.round( delta_nonlin_wk.data[:], decimals=6) )
    assert np.shape(delta_nonlin) == (1,)
    delta_nonlin = np.abs(delta_nonlin[0])
    
    # Compare with gap function evaluated from analytical Matsubara summations
    anaGapEqn = AnalyticalMatsubaraSumGapEquation(delta_nonlin, beta, q_fmp, U, mu)
    if leadingEv >= 1.0:
        np.testing.assert_array_almost_equal(anaGapEqn, 1.0, decimal=3)
    
    return leadingEv, delta_nonlin

def test_nonlinear_Eliashberg_solver(verbose=False):
    """ Test the non-linear Eliahsberg solver against an analytical result and the linearized Eliashberg solver
    Author: Yann in 't Veld (2024) """

    fmp_index = 1
    
    if verbose:
        TLst = [200, 300, 400, 500, 600]
    else:
        TLst = [200, 600]
    datLst = np.array([SolveEliashbergAtTemperature(T, fmp_index) for T in TLst])

    leadingEvLst = datLst[:,0]
    delta_nonlin_Lst = datLst[:,1]
            
    if verbose:
        import matplotlib.pyplot as plt

        leadingEvLst = np.array(leadingEvLst, dtype=float)
        Tc = np.interp(1.0, leadingEvLst[::-1], TLst[::-1])
        print("Tc =", Tc, "K")

        fig = plt.figure(figsize=(5.0,3.0),dpi=300)

        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()

        ax1.plot(TLst, leadingEvLst, marker='.', color='cornflowerblue')
        ax1.axhline(y=1.0, linewidth=0.5, color='gray', linestyle=':')
        ax1.set_ylabel("$\lambda$", color="cornflowerblue")

        ax2.plot(TLst, delta_nonlin_Lst, marker='.', color='firebrick')
        ax2.axhline(y=0.0, linewidth=0.5, color='black', linestyle='--')
        ax2.set_ylabel("$\Delta$", color="firebrick")

        ax1.set_xlabel("T (K)")

        plt.tight_layout()
        plt.show()
        
if __name__ == "__main__":
    test_nonlinear_Eliashberg_solver(verbose=False)
