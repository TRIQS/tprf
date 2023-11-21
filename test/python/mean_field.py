""" Mean field and ferro magnetic q=0 response for the one dimensional 
Hubbard model. 

Comparing numerical calculation with rank 4 tensor and analytical 
results from Fazekas (Ch. 7.4)

Author: Hugo U. R. Strand (2018) """

# ----------------------------------------------------------------------

import numpy as np

from scipy.integrate import quad
from scipy.optimize import brentq
    
# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq, Idx

# ----------------------------------------------------------------------

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lindhard_chi00

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------
def get_kinetic_energy_ref(t, beta):

    def disp(k):
        return -2*t * np.cos(k)
    
    def integrand(k):
        e = disp(k)
        f = 1./( np.exp(beta * e) + 1 )
        return e * f

    E_kin, err = quad(integrand, -np.pi, np.pi)

    n_spin = 2.
    k_vol = 2*np.pi
    
    E_kin *= n_spin / k_vol

    return E_kin

# ----------------------------------------------------------------------
def get_total_energy_mf_ref(t, beta, U, mu, n, m):
    E_kin = get_kinetic_energy_mf_ref(t, beta, U, mu, n, m)
    E_tot = E_kin - U*(0.25*n**2 - m**2)
    return E_tot
    
# ----------------------------------------------------------------------
def get_kinetic_energy_mf_ref(t, beta, U, mu, n, m):

    def disp(k, U, n, m, s):
        return -2*t * np.cos(k) - s*U*m + 0.5*U*n
    
    def integrand(k, U, mu, n, m, s):
        e = disp(k, U, n, m, s)
        f = 1./( np.exp(beta * (e - mu)) + 1 )
        return e * f

    E_kin_up, err = quad(integrand, -np.pi, np.pi, args=(U, mu, n, m, +1.))
    E_kin_do, err = quad(integrand, -np.pi, np.pi, args=(U, mu, n, m, -1.))

    E_kin = E_kin_up + E_kin_do
    
    k_vol = 2*np.pi    
    E_kin *= 1. / k_vol

    return E_kin

# ----------------------------------------------------------------------
def get_kinetic_energy(e_k, g_wk):

    E_kin = 0.

    kmesh = e_k.mesh
    
    for k in kmesh:
        g_w = g_wk[:, k][0,0].copy()
        g_w.data[:] = np.einsum('wab,ab->w', g_wk[:, k].data, e_k[k])
        E_kin += g_w.density()

    E_kin /= len(e_k.mesh)

    return E_kin.real

# ----------------------------------------------------------------------
def get_density_matrix(g_wk):

    kmesh = g_wk.mesh.components[1]
    
    rho = 0.
    for k in kmesh:
        rho += g_wk[:, k].density()

    rho /= len(kmesh)

    return rho.real

# ----------------------------------------------------------------------
def E_vs_m(t, beta, U):

    mu = 0.5*U

    n = 1.0

    m_vec = np.linspace(-0.5, 0.5, num=100)

    E_vec = np.array(
        [ get_total_energy_mf_ref(t, beta, U, mu, n, m) for m in m_vec ])

    return m_vec, E_vec

# ----------------------------------------------------------------------
def fermi_distribution(beta, eps):
    return 1./(np.exp(beta*eps) + 1.)

# ----------------------------------------------------------------------
def fermi_distribution_derivative(beta, eps):
    return -beta/4. / np.cosh(0.5*beta*eps)**2

# ----------------------------------------------------------------------
def get_tb_model(t, U, n, m, mu=0., vanilla_tb=False):
    
    h_loc = -U*m*np.diag([+1., -1.]) + (0.5*U*n - mu) * np.eye(2)
    T = - t * np.eye(2)

    if vanilla_tb:
        from triqs.lattice.tight_binding import TBLattice
    else:
        from triqs_tprf.tight_binding import TBLattice
        
    t_r = TBLattice(
        units = [(1, 0, 0)],
        hopping = {
            # nearest neighbour hopping -t
            (0,): h_loc,
            (+1,): T,
            (-1,): T,
            },
        orbital_positions = [(0,0,0)] * 2,
        orbital_names = ['up', 'do'],
        )

    return t_r

# ----------------------------------------------------------------------
def test_fermi_distribution_derivative(verbose=False):

    beta = 3.3
    e = np.linspace(-5., 5., num=1000)
    f = fermi_distribution(beta, e)
    df = fermi_distribution_derivative(beta, e)

    from scipy.interpolate import InterpolatedUnivariateSpline
    sp = InterpolatedUnivariateSpline(e, f)
    dsp = sp.derivative(1)
    df_ref = dsp(e)

    np.testing.assert_almost_equal(df, df_ref)

    if verbose:
        import matplotlib.pyplot as plt
        plt.plot(e, f, label='f')
        plt.plot(e, df, label='df')
        plt.plot(e, df_ref, label='df_ref')
        plt.legend()
        plt.show()
        exit()

# ----------------------------------------------------------------------
def get_density_of_states(eps, t):

    return 1./np.sqrt(1. - (eps/(2*t))**2) / (2*t * np.pi)
    
# ----------------------------------------------------------------------
def test_tb_model(verbose=False):

    t = 1.0
    U = 0.0
    mu = 0.5*U

    m = 0.0
    n = 1.0
    
    t_r = get_tb_model(t, U, n, m, mu=0., vanilla_tb=True)

    from triqs.lattice.tight_binding import dos
    rhos = dos(t_r.tb, n_kpts=100000, n_eps=100, name='foobar')

    eps = rhos[0].eps
    rho_ref = get_density_of_states(eps, t)

    np.testing.assert_almost_equal(rhos[0].rho[1:-1], rho_ref[1:-1], decimal=2)

    if verbose:
        import matplotlib.pyplot as plt    
        for rho in rhos:
            plt.plot(rho.eps, rho.rho)

        plt.plot(eps, rho_ref)
        plt.ylim([0, 4])
        plt.show()
        exit()
        
# ----------------------------------------------------------------------
def chi0_q0_integral(t, beta):

    def integrand(eps):
        rho = get_density_of_states(eps, t)
        df = fermi_distribution_derivative(beta, eps)
        return -rho * df

    chi0_q0, err = quad(integrand, -2.*t, 2.*t)

    return chi0_q0        

# ----------------------------------------------------------------------
def find_Uc(t, beta):
    
    def root_function(U):
        chi0_q0 = chi0_q0_integral(t, beta)
        return 1 - U * chi0_q0

    Uc = brentq(root_function, 0, 100.)

    print('Uc =', Uc)

    return Uc

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    test_fermi_distribution_derivative(verbose=False)
    test_tb_model(verbose=False)
    
    n_k = (128, 1, 1)
    nw = 400

    beta = 5.0

    t = 1.0
    U = 6.1
    mu = 0.5*U

    m = 0.0
    n = 1.0
    
    E_kin_ref = get_kinetic_energy_mf_ref(t, beta, U, mu, n, m)
    print('E_kin_ref =', E_kin_ref)

    E_tot_ref = get_total_energy_mf_ref(t, beta, U, mu, n, m)
    print('E_tot_ref =', E_tot_ref)

    print('--> tight binding model')
    t_r = get_tb_model(t, U, n, m, mu=0.)

    print('--> dispersion e_k')
    kmesh = t_r.get_kmesh(n_k)
    e_k = t_r.fourier(kmesh)

    #print e_k.data
    
    print('--> lattice g0_wk')
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

    E_kin = get_kinetic_energy(e_k, g0_wk)
    print('E_kin =', E_kin)
    
    np.testing.assert_almost_equal(E_kin_ref, E_kin, decimal=6)

    rho = get_density_matrix(g0_wk)
    print('rho =\n', rho)

    n = rho[0, 0] + rho[1, 1]
    m = 0.5 * (rho[0, 0] - rho[1, 1])

    print('n, m =', n, m)

    E_tot = E_kin - U*(n**2/4 - m**2)
    print('E_tot =', E_tot)

    np.testing.assert_almost_equal(E_tot_ref, E_tot, decimal=6)

    # ------------------------------------------------------------------
    # -- Lattice chi0
    
    print('--> chi00_wk')
    chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)
    print('chi0_q0 =\n', chi00_wk[Idx(0), Idx(0, 0, 0)].real.reshape((4,4)))

    print('--> lindhard_chi00')
    wmesh = MeshImFreq(beta, "Boson", 1)
    chi00_wk_analytic = lindhard_chi00(e_k=e_k, mesh=wmesh, mu=mu)
    print('chi0_q0_analytic =\n', chi00_wk_analytic[
        Idx(0), Idx(0, 0, 0)].real.reshape((4,4)))

    np.testing.assert_almost_equal(
        chi00_wk.data, chi00_wk_analytic.data, decimal=5)

    chi0_q0_ref = chi0_q0_integral(t, beta)
    
    print('chi0_q0     =', chi00_wk[Idx(0), Idx(0, 0, 0)][0,0,0,0].real)
    print('chi0_q0_ref =', chi0_q0_ref)

    np.testing.assert_almost_equal(
        chi00_wk_analytic[Idx(0), Idx(0, 0, 0)][0,0,0,0], chi0_q0_ref)

    # ------------------------------------------------------------------
    # -- RPA tensor
    
    gf_struct = [[0, 2]]
    
    from triqs.operators import n, c, c_dag, Operator, dagger
    H_int = U * n(0, 0) * n(0, 1)
    print('H_int =', H_int)

    from triqs_tprf.rpa_tensor import get_rpa_tensor
    from triqs_tprf.rpa_tensor import fundamental_operators_from_gf_struct

    fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)
    print(fundamental_operators)
    U_abcd = get_rpa_tensor(H_int, fundamental_operators)

    print('U_abcd =\n', U_abcd.reshape((4, 4)))
    
    # ------------------------------------------------------------------
    # -- RPA

    from triqs_tprf.lattice import solve_rpa_PH

    chi_wk = solve_rpa_PH(chi00_wk, U_abcd)
    print('chi_q0 =\n', chi_wk[Idx(0), Idx(0, 0, 0)].real.reshape((4,4)))

    Sz = 0.5 * np.diag([+1., -1.])
    chi_SzSz_wk = chi_wk[0,0,0,0].copy()
    chi_SzSz_wk.data[:] = np.einsum('wkabcd,ab,cd->wk', chi_wk.data, Sz, Sz)
    print('chi_SzSz_q0 =', chi_SzSz_wk[Idx(0), Idx(0, 0, 0)].real)

    # Eq. 7.43 Fazekas (additional 0.5 factor)
    chi_q0_ref = 0.5 * chi0_q0_ref / (1. - U*chi0_q0_ref)
    print('chi_q0_ref =', chi_q0_ref)

    np.testing.assert_almost_equal(
        chi_SzSz_wk[Idx(0), Idx(0, 0, 0)], chi_q0_ref)
    
    # ------------------------------------------------------------------
    # -- Mean field results for E[m] and Uc
    
    if False:
        
        import matplotlib.pyplot as plt
        plt.figure(figsize=(3.25*2, 7))
        
        t = 1.0
        beta = 5.0
        U_vec = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8.])
        
        for U in U_vec:
            
            m, E = E_vs_m(t, beta, U)
            chi0_q0 = chi0_q0_integral(t, beta)

            subp = [2, 1, 1]
            
            plt.subplot(*subp); subp[-1] += 1
            plt.plot(m, E, '-', lw=0.5, label='U=%2.2f' % U)

            plt.subplot(*subp); subp[-1] += 1
            plt.plot(U, 1 - chi0_q0 * U, 'x')
            plt.grid()
            
        Uc = find_Uc(t, beta)
        m, E = E_vs_m(t, beta, Uc)
        subp = [2, 1, 1]
            
        plt.subplot(*subp); subp[-1] += 1
        plt.plot(m, E, '-k', label='U=%2.2f' % Uc)
        plt.xlabel(r'$m$')
        plt.ylabel(r'$E[m]$')
        plt.legend()
        
        plt.subplot(*subp); subp[-1] += 1
        plt.plot(Uc, 0., '+r', label='Uc=%2.2f' % Uc)
        plt.ylabel(r'$1 - U \chi_0(q=0)$')
        plt.xlabel(r'$U$')
        plt.legend()

        plt.tight_layout()
        plt.savefig('figure_E_vs_m_and_Uc.pdf')
        plt.show()

        exit()
