# ----------------------------------------------------------------------

import numpy as np

from triqs_tprf.lattice import lattice_dyson_g0_wk
#from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.gw import gw_sigma, g0w_sigma

from triqs.gf import Gf, MeshImFreq, MeshReFreq, MeshBrZone
from triqs.gf.mesh_product import MeshProduct
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

def plot_linecut_w(ax, g_fk, kpoint, target=(0,0), opt={}, wmin=-np.inf, wmax=np.inf, adj=lambda x:x):
    fmesh, kmesh = g_fk.mesh.components
 
    A0_fk_interper = np.vectorize(
        lambda w, kx, ky, kz : adj(g_fk[w,:]([kx, ky, kz])[target]))
    
    A0_fk_interp = []
    farr = []
    for w in fmesh:
        if(w.value < wmin or w.value > wmax): continue
        farr += [w.value]
        A0_fk_interp += [A0_fk_interper(w, kpoint[0], kpoint[1], kpoint[2])]
    
    p = ax.plot(farr, A0_fk_interp, **opt)
    return p

def nF(ww, beta):
    """Fermi-Dirac distribution function"""
    return 0.5 - 0.5 * np.tanh(0.5 * beta * ww)

def nB(ww, beta):
    """Bose-Einstein distribution function"""
    return 1.0 / np.expm1( beta * ww )

def ExactSigma0D(iw, beta, g2, wD, E):
    """Exact result for the 0-dimensional GW self-energy for an electron-phonon interaction with dispersionless
    dispersion and scalar electron-phonon coupling, using a single k-point at gamma."""
    fact1 = (nB(wD, beta) + nF(E, beta)) / (iw + wD - E)
    fact2 = (nB(wD, beta) + 1.0 - nF(E, beta)) / (iw - wD - E)
    return g2 * (fact1 + fact2)

def ExactSigma0D_multiband(iw, beta, g2, wD, Eorb):
    """Exact result for the 0-dimensional GW self-energy for an electron-phonon interaction with dispersionless
    dispersion and scalar electron-phonon coupling, using a single k-point at gamma. For a multi-band system."""

    Eband, T = np.linalg.eig(Eorb)
    assert np.allclose(T@np.diag(Eband)@np.linalg.inv(T), Eorb)

    norb = len(Eband)
    sigma = np.zeros((norb, norb), dtype=complex)
    for a in range(norb):
        for b in range(norb):
            for l in range(norb):
                sigma[a,b] += T[a,l] * np.linalg.inv(T)[l,b] * ExactSigma0D(iw, beta, g2[a,b], wD[a,b], Eband[l])
    return sigma
    
    
def ElectronPhononInteraction(iw, g2, wD):
    """Electron-phonon interaction with a dispersionless phonon wD and a scalar electron-phonon coupling g2"""
    return g2 * 2.0 * wD / (iw**2.0 - wD**2.0)

def test_gw_sigma_against_exact_Matsubara():
    """ Tests the Matsubara frequency axis GW implementation by comparing to an exact analytic result.
    This result was found for a calculation on a single k-point, with a simple electron-phonon propagator.
    Author: Yann in 't Veld (2023) """ 
    
    print("== Matsubara axis ==")

    g2 = np.array([[5, 0.5],
                   [0.5, 5]])
    wD = np.array([[0.5, 1],
                   [1, 1.5]])
    Eorb = np.array([[5, 2],
                     [2,-3]])
    beta = 300.0
    nw = 500
    wmin = -10
    wmax = 10
    eta = 1e-2
    norb = 2
    nw = 500

    print('--> construct mesh and Enk')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)]*norb)
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, np.array([1, 1, 1], dtype=int))
    fmesh = MeshReFreq(wmin, wmax, nw)
    
    Enk = Gf(mesh=kmesh, target_shape=[norb]*2)   
    Enk.data[0,:,:] = Eorb
    
    wmesh = MeshImFreq(beta, 'Fermion', nw)
    numesh = MeshImFreq(beta, 'Boson', nw)

    print('--> lattice_dyson_g0_wk')
    g0_fk = lattice_dyson_g0_wk(0.0, Enk, wmesh)

    print('--> bare electron-phonon interaction')
    I_phon_wk = Gf(mesh=MeshProduct(numesh, kmesh), target_shape=[norb]*4)
    I_phon_wk.data[:] = 0.0
    for nu in numesh:
        nuii = nu.data_index
        for a in range(norb):
            for b in range(norb):
                I_phon_wk.data[nuii,0,a,a,b,b] =  ElectronPhononInteraction(nu.value, g2[a,b], wD[a,b])

    I_phon_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    I_phon_k.data[:] = 0.0 

    print('--> gw_sigma')
    sigma_wk = gw_sigma(I_phon_wk, g0_fk)

    sigma_ref_wk = Gf(mesh=sigma_wk.mesh, target_shape=sigma_wk.target_shape)
    for f in wmesh:
        fii = f.data_index
        sigma_ref_wk.data[fii,:] = ExactSigma0D_multiband(f.value, beta, g2, wD, Eorb)

    np.testing.assert_array_almost_equal(sigma_wk.data[:], sigma_ref_wk.data[:], decimal=1e-6)


def test_gw_sigma_against_exact_realfreq(verbose=False):
    """ Tests thereal-frequency axis GW implementation by comparing to an exact analytic result.
    This result was found for a calculation on a single k-point, with a simple electron-phonon propagator.
    Author: Yann in 't Veld (2023) """ 

    print("== Real-freq axis ==")

    g2 = np.array([[5, 0.5],
                   [0.5, 5]])
    wD = np.array([[0.5, 1],
                   [1, 1.5]])
    Eorb = np.array([[5, 2],
                     [2,-3]])
    beta = 300.0
    nw = 500
    wmin = -10
    wmax = 10
    eta = 1e-2
    norb = 2

    print('--> construct mesh and Enk')
    bl = BravaisLattice(units=[(1,0,0)], orbital_positions=[(0,0,0)]*norb)
    bz = BrillouinZone(bl)
    kmesh = MeshBrZone(bz, np.array([1, 1, 1], dtype=int))
    fmesh = MeshReFreq(wmin, wmax, nw)
    
    Enk = Gf(mesh=kmesh, target_shape=[norb]*2)   
    Enk.data[0,:,:] = Eorb
    
    
    fmesh = MeshReFreq(wmin, wmax, nw)

    print('--> bare electron-phonon interaction')
    I_phon_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=[norb]*4)
    for f in fmesh:
        fii = f.data_index
        for a in range(norb):
            for b in range(norb):
                I_phon_fk.data[fii,0,a,a,b,b] =  ElectronPhononInteraction(f.value + 1.0j*eta, g2[a,b], wD[a,b])

    I_phon_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    I_phon_k.data[:] = 0.0 

    if verbose:
        import matplotlib.pyplot as plt

        print('--> plotting I_fk')
        fig = plt.figure(figsize=(12.0,3.0),dpi=300)
        
        ax = fig.add_subplot(141)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x= wD[0,0], color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=-wD[0,0], color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(0,0,0,0), opt={"label":r"Re$(I)$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(0,0,0,0), opt={"label":r"Im$(I)$", "linewidth":1}, adj=lambda x:np.imag(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.set_ylabel(r"$I_{k=\Gamma}(\omega)$ (eV)")
        ax.legend()
        ax.set_title("Element = [0,0]")
        
        ax = fig.add_subplot(142)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x= wD[1,1], color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=-wD[1,1], color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(1,1,1,1), opt={"label":r"Re$(I)$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(1,1,1,1), opt={"label":r"Im$(I)$", "linewidth":1}, adj=lambda x:np.imag(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.legend()
        ax.set_title("Element = [1,1]")
        
        ax = fig.add_subplot(143)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x= wD[0,1], color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=-wD[0,1], color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(0,0,1,1), opt={"label":r"Re$(I)$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(0,0,1,1), opt={"label":r"Im$(I)$", "linewidth":1}, adj=lambda x:np.imag(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.legend()
        ax.set_title("Element = [0,1]")
        
        ax = fig.add_subplot(144)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x= wD[1,0], color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=-wD[1,0], color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(1,1,0,0), opt={"label":r"Re$(I)$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, I_phon_fk, [0,0,0], target=(1,1,0,0), opt={"label":r"Im$(I)$", "linewidth":1}, adj=lambda x:np.imag(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.legend()
        ax.set_title("Element = [1,0]")
        
        plt.show()
    
    print('--> gw_sigma')
    sigma_fk = g0w_sigma(0.0, beta, Enk, I_phon_fk, I_phon_k, eta)

    print('--> reference sigma')
    sigma_ref_fk = Gf(mesh=sigma_fk.mesh, target_shape=sigma_fk.target_shape)
    for f in fmesh:
        fii = f.data_index
        sigma_ref_fk.data[fii,0,:] = ExactSigma0D_multiband(f.value + 1.0j*eta, beta, g2, wD, Eorb)

    if verbose:
        print('--> plotting sigma_fk')
        Eband, _ = np.linalg.eig(Eorb)

        fig = plt.figure(figsize=(12.0,3.0),dpi=300)
        
        ax = fig.add_subplot(141)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x=Eband[0].real, color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=Eband[1].real, color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, sigma_fk, [0,0,0], target=(0,0), opt={"label":r"Re$(\Sigma^{num})$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, sigma_ref_fk, [0,0,0], target=(0,0), opt={"label":r"Re$(\Sigma^{ana})$", "linewidth":1, "linestyle":"--"}, adj=lambda x:np.real(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.set_ylabel(r"$\Sigma_{k=\Gamma}(\omega)$ (eV)")
        ax.legend()
        ax.set_title("Element [0,0]")
        
        ax = fig.add_subplot(142)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x=Eband[0].real, color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=Eband[1].real, color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, sigma_fk, [0,0,0], target=(1,1), opt={"label":r"Re$(\Sigma^{num})$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, sigma_ref_fk, [0,0,0], target=(1,1), opt={"label":r"Re$(\Sigma^{ana})$", "linewidth":1, "linestyle":"--"}, adj=lambda x:np.real(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.legend()
        ax.set_title("Element [1,1]")
        
        ax = fig.add_subplot(143)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x=Eband[0].real, color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=Eband[1].real, color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, sigma_fk, [0,0,0], target=(0,1), opt={"label":r"Re$(\Sigma^{num})$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, sigma_ref_fk, [0,0,0], target=(0,1), opt={"label":r"Re$(\Sigma^{ana})$", "linewidth":1, "linestyle":"--"}, adj=lambda x:np.real(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.legend()
        ax.set_title("Element [0,1]")
        
        ax = fig.add_subplot(144)
        ax.axhline(y=0.0, color='black',linewidth=0.5,linestyle='-')
        ax.axvline(x=Eband[0].real, color='black', linewidth=0.5,linestyle=':')
        ax.axvline(x=Eband[1].real, color='black', linewidth=0.5,linestyle=':')
        plot_linecut_w(ax, sigma_fk, [0,0,0], target=(1,0), opt={"label":r"Re$(\Sigma^{num})$", "linewidth":1}, adj=lambda x:np.real(x))
        plot_linecut_w(ax, sigma_ref_fk, [0,0,0], target=(1,0), opt={"label":r"Re$(\Sigma^{ana})$", "linewidth":1, "linestyle":"--"}, adj=lambda x:np.real(x))
        ax.set_xlabel(r"$\omega$ (eV)")
        ax.legend()
        ax.set_title("Element [1,0]")
        
        plt.show()     

    print('--> compare')
    # Rather large tolerance needed. Probably due to the numerical spectral function of I not being a
    # true delta-function, and thus deviating from the analytic result.
    assert np.allclose(sigma_ref_fk.data[:], sigma_fk.data[:], rtol=10.0)
    
    
if __name__ == "__main__":
    test_gw_sigma_against_exact_Matsubara()
    test_gw_sigma_against_exact_realfreq(verbose=False)
