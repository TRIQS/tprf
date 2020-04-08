# ----------------------------------------------------------------------

""" One dimensional Hubbard model solved with Hartree-Fock

Compute the FM spin susceptibility by RPA, -dm/dh, and -d^2Omega/dh^2
and compare the results.

Author: Hugo U.R. Strand (2018), hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline

# ----------------------------------------------------------------------

from pytriqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.hf_solver import HartreeSolver
from triqs_tprf.hf_response import HartreeResponse

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    beta = 1.3
    N_tot = 0.7
    n_k = (256, 1, 1)

    # -- One dimensional tight binding model

    t = 1.0
    h_loc = np.zeros((2, 2))
    T = - t * np.eye(2)
        
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

    e_k = t_r.on_mesh_brillouin_zone(n_k)

    # -- Local double occupancy operator
    
    gf_struct = [[0, [0, 1]]]

    docc = n(0, 0) * n(0, 1)

    Sz = 0.5 * np.diag([1., -1.])

    # -- Sweep in interaction U
    
    U_vec = np.arange(9., -0.5, -1.5)
    h_vec = 1e-3 * np.linspace(-1., 1., num=11)
    
    res = []

    M0 = np.zeros((2, 2))
    M = -5.0 * Sz
    mu = 0.
    
    for U in U_vec:

        print('-'*72)
        print('U =', U)
        print('-'*72)

        mu0 = 0.5 * U # half-filling
        
        H_int = U * docc

        hs = HartreeSolver(e_k, beta, H_int=H_int, gf_struct=gf_struct)

        # -- Solve without seeding FM symmetry breaking
        hs.solve_newton(N_target=N_tot, M0=M0, mu0=mu0)

        # -- Compute the RPA response of the non-magnetic solution
        hr = HartreeResponse(hs)
        hs.chi = hr.response(Sz, Sz)

        # -- Apply magnetic field

        hs.m_vec = []
        hs.E_vec = []
        hs.E_kin_vec = []
        hs.E_int_vec = []
        hs.Omega0_vec = []
        hs.Omega_vec = []
        hs.h_vec = h_vec
        
        for h in hs.h_vec:
            e_k_mag = e_k.copy()
            e_k_mag.data[:] += h * Sz[None, ...]

            hsh = HartreeSolver(e_k_mag, beta, H_int=H_int, gf_struct=gf_struct)
            hsh.solve_newton(N_target=N_tot, M0=hs.M, mu0=hs.mu)

            hsh.m = hsh.expectation_value(Sz)
            hs.m_vec.append(hsh.m)

            hs.E_vec.append(hsh.E_tot)
            hs.E_kin_vec.append(hsh.E_kin)
            hs.E_int_vec.append(hsh.E_int)

        for h in hs.h_vec:
            e_k_mag = e_k.copy()
            e_k_mag.data[:] += h * Sz[None, ...]

            hsh = HartreeSolver(e_k_mag, beta, H_int=H_int, gf_struct=gf_struct)
            hsh.solve_newton_mu(hs.mu, M0=hs.M)
            
            hs.Omega0_vec.append(hsh.Omega0)
            hs.Omega_vec.append(hsh.Omega)
            
        # -- Compute susceptibility from field response

        spl_m_h = InterpolatedUnivariateSpline(hs.h_vec, hs.m_vec)
        hs.chi_dmdh = -spl_m_h(0., nu=1)

        spl_Omega_h = InterpolatedUnivariateSpline(hs.h_vec, hs.Omega_vec)
        hs.chi_d2Omegadh2 = -spl_Omega_h(0., nu=2)

        print('chi            =', hs.chi)
        print('chi_dmdh       =', hs.chi_dmdh)
        print('chi_d2Omegadh2 =', hs.chi_d2Omegadh2)

        #exit()
                
        # -- Solve by seeding FM symmetry breaking
        
        hs.solve_newton(N_target=N_tot, M0=M, mu0=mu)
        mu = hs.mu
        M = hs.M.copy()

        hs.m = hs.expectation_value(Sz)
        hs.U = U
        
        res.append(hs)            
        
    get_vec = lambda vec, attr : np.array([ getattr(el, attr) for el in vec])

    U = get_vec(res, 'U')
    m = get_vec(res, 'm')
    chi = get_vec(res, 'chi')
    chi_dmdh = get_vec(res, 'chi_dmdh')
    chi_d2Omegadh2 = get_vec(res, 'chi_d2Omegadh2')

    E_tot = get_vec(res, 'E_tot')
    E_kin = get_vec(res, 'E_kin')
    E_int = get_vec(res, 'E_int')

    Omega0 = get_vec(res, 'Omega0')
    Omega = get_vec(res, 'Omega')
    
    m_vec = get_vec(res, 'm_vec')
    E_vec = get_vec(res, 'E_vec')
    E_kin_vec = get_vec(res, 'E_int_vec')
    E_int_vec = get_vec(res, 'E_kin_vec')
    Omega0_vec = get_vec(res, 'Omega0_vec')
    Omega_vec = get_vec(res, 'Omega_vec')
    h_vec = res[0].h_vec

    chi_diff1 = 1./chi - 1./chi_dmdh
    chi_diff2 = 1./chi - 1./chi_d2Omegadh2

    print('chi_diff1 =', np.max(np.abs(chi_diff1)))
    print('chi_diff2 =', np.max(np.abs(chi_diff2)))

    np.testing.assert_array_almost_equal(1./chi, 1./chi_dmdh)
    np.testing.assert_array_almost_equal(1./chi, 1./chi_d2Omegadh2, decimal=3)

    if False:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(7, 7))

        subp = [3, 2, 1]
        plt.subplot(*subp); subp[-1] += 1

        plt.plot(U, m, '.-')
        plt.xlabel(r'$U$')
        plt.ylabel(r'$m$')
        plt.subplot(*subp); subp[-1] += 1

        plt.plot(U, 1./chi, 'o-', label=r'HF-RPA $[\chi]^{-1}$')
        plt.plot(U, 1./chi_dmdh, 'x-',
                 label=r'$[-\frac{\partial m}{\partial h}|_{\langle \hat{N} \rangle = N, h=0}]^{-1}$')
        plt.plot(U, 1./chi_d2Omegadh2, '.-',
                 label=r'$[-\frac{\partial^2\Omega}{\partial h^2}|_{\mu, h=0}]^{-1}$')

        plt.xlabel(r'$U$')
        plt.ylabel(r'$1/\chi$')
        plt.legend(loc='best')
        plt.grid()

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(U, E_tot, label='E_tot')
        plt.plot(U, E_kin, label='E_kin')
        plt.plot(U, E_int, label='E_int')
        plt.legend(loc='best')
        plt.ylabel(r'$U$')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(U, Omega0, label=r'$\Omega_0$')
        plt.plot(U, Omega, label=r'$\Omega$')
        plt.legend(loc='best')
        plt.ylabel(r'$U$')

        plt.subplot(*subp); subp[-1] += 1
        for u, m in zip(U, m_vec):
            plt.plot(h_vec, m, '-', label='U = %2.2f' % u)

        #plt.legend(loc='best')
        plt.xlabel(r'$h$')
        plt.ylabel(r'$m$')

        plt.subplot(*subp); subp[-1] += 1

        for u, O in zip(U, Omega_vec):
            O_h0 = InterpolatedUnivariateSpline(h_vec, O)(0.)
            plt.plot(h_vec, O - O_h0, '-')

        plt.xlabel(r'$h$')
        plt.ylabel(r'$\Omega$')

        plt.tight_layout()
        plt.savefig('figure_1d_hubbard_hf_rpa.pdf')
        plt.show()

