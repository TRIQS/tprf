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

def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)

def InitialConditionDependenceTest(verbose):     
    """ Test of the non-linear Eliashberg solver gives the same result fordifferent initial guesses.
    Author: Yann in 't Veld (2024) """

    mu = 1.0
    g2 = 0.3
    wD = 0.3
    nk = 8
    U = 0.5
    T = 200.0
    wmax = 50.
    eps = 1e-10
    
    
    print('--> construct meshes')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)])
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, [nk, 1, 1])
    
    beta = 1.0 / (kB * T)
    wmesh_f  = MeshDLRImFreq(beta, 'Fermion', wmax, eps)
    wmesh_b = MeshDLRImFreq(beta, 'Boson', wmax, eps)
    
    Enk = Gf(mesh=kmesh, target_shape=[1]*2)
    for k in kmesh:
        kfolded = backfold_k(k.value, kmesh)
        knorm = np.linalg.norm(kfolded)
        Enk.data[k.data_index,:] = 0.1 * knorm**2.0

    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=Enk, mesh=wmesh_f)

    print('--> setup interaction vertex')
    I_k = Gf(mesh=kmesh, target_shape=[1]*4)
    I_k.data[:] = 0.0
    for k in kmesh:
        I_k.data[:] = U

    I_wk = Gf(mesh=MeshProduct(wmesh_b, kmesh), target_shape=[1]*4)
    for w in wmesh_b:
        wii = w.data_index
        I_wk.data[wii,:] = ElectronPhononInteraction(w.value, g2 ,wD) + I_k.data[:]
        
        
    print("--> solve_eliashberg_nonlinear (start from eigenvector of linearized Eliashberg)")
    delta0_wk = semi_random_initial_delta(g0_wk)
    vals, vecs = solve_eliashberg(I_wk, g0_wk, initial_delta=delta0_wk, Gamma_pp_const_k=I_k, product="FFT", solver="IRAM", tol=1e-8)
    leadingIndex = np.argmax(np.real(vals))
    leadingEv = vals[leadingIndex]
    delta_lin_wk = vecs[leadingIndex]
    print("    leading ev:", leadingEv)
    
    delta_nonlin_wk = solve_eliashberg_nonlinear(I_wk, g0_wk, initial_delta=delta_lin_wk, Gamma_pp_const_k=I_k, product="FFT", tol=1e-8)
    delta_nonlin_wk.data[:] = np.abs(delta_nonlin_wk.data[:])
        
    
    
    print("--> solve_eliashberg_nonlinear (start from randomly generated intial guesses)")

    seeds = np.random.randint(0, 1000, size=10)
    deltaLst = []
    for seed in seeds:
        print("Seed:", seed)
        delta0_wk = semi_random_initial_delta(g0_wk, seed=seed)
        delta_wk = solve_eliashberg_nonlinear(I_wk, g0_wk, initial_delta=delta0_wk, Gamma_pp_const_k=I_k, product="FFT", tol=1e-8)
        delta_wk.data[:] = np.abs(delta_wk.data[:])
        deltaLst += [delta_wk]
        
        np.testing.assert_array_almost_equal(np.abs(delta_wk.data[:]), np.abs(delta_nonlin_wk.data[:]))
        
        
    if verbose:
        from triqs.plot.mpl_interface import oplot,plt
        
        fig = plt.figure(figsize=(5.0,3.0),dpi=300)
        ax = fig.add_subplot(111)

        delta_nonlin_w = Gf(mesh=wmesh_f, target_shape=delta_nonlin_wk.target_shape)
        delta_nonlin_w.data[:] = np.sum(delta_nonlin_wk.data[:], axis=1) / len(kmesh)
        oplot(delta_nonlin_w.real, linewidth=1.0, linestyle=':', label="Lin. $\Delta_0$")

        for ii,delta_wk in enumerate(deltaLst):
            delta_w = Gf(mesh=wmesh_f, target_shape=delta_wk.target_shape)
            delta_w.data[:] = np.sum(delta_wk.data[:], axis=1) / len(kmesh)
            oplot(delta_w.real, linewidth=1.0, linestyle='-', marker='.', label=ii)

        ax.set_xlim(-4, 4)
        ax.set_xlabel("$i\omega_m$")
        ax.set_ylabel("$|\Delta(i\omega_{n})|$")
        plt.show()
    
if __name__ == "__main__":
    InitialConditionDependenceTest(verbose=False)
