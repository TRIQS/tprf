# ----------------------------------------------------------------------

""" One dimensional Hubbard model solved with Hartree-Fock

Compute the AF and FM solutions for a two-site super cell
and compare the transition points to the divergencies of
the respective, susceptibility.

Author: Hugo U.R. Strand (2018), hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline

# ----------------------------------------------------------------------

from triqs.gf import Idx
from triqs.lattice.lattice_tools import BrillouinZone
from triqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.super_lattice import TBSuperLattice

from triqs_tprf.hf_solver import HartreeFockSolver
from triqs_tprf.hf_response import HartreeFockResponse

# ----------------------------------------------------------------------
if __name__ == '__main__':

    np.set_printoptions(
        precision=3, linewidth=10000,
        threshold=1000000, suppress=True)

    beta = 1.2
    N_tot = 0.7 * 2
    n_k = (8, 1, 1)

    # -- One dimensional tight binding model

    t = 1.0
    h_loc = np.zeros((2, 2))
    T = - t * np.eye(2)
        
    t_r_prim = TBLattice(
        units = [(1,)],
        hopping = {
            # nearest neighbour hopping -t
            (0,): h_loc,
            (+1,): T,
            (-1,): T,
            },
        orbital_positions = [(0,0,0)] * 2,
        orbital_names = ['up', 'do'],
        )    

    # -- Two site super lattice
    
    P = np.array([[ 2 ]]) 
    
    units_prim = np.array(t_r_prim.units)
    print('units_prim =\n', units_prim)
    units = np.dot(P, units_prim)
    print('units =\n', units)
    
    t_r = TBSuperLattice(t_r_prim, P)

    kmesh = t_r.get_kmesh(n_k)
    e_k = t_r.fourier(kmesh)
    print(e_k.target_shape)
    
    print('eps(k=0) =\n', e_k[Idx(0, 0, 0)])

    # -- Local double occ and spin operators

    gf_struct = [[0, 4]]

    docc = n(0, 0) * n(0, 2) + n(0, 1) * n(0, 3)

    print('docc =', docc)
    
    sigma_x = 0.5 * np.rot90(np.diag([1., 1.]))
    sigma_y = 0.5 * np.rot90(np.diag([1.j, -1.j]))
    sigma_z = 0.5 * np.diag([1., -1.])

    print('sigma_x =\n', sigma_x)
    print('sigma_y =\n', sigma_y)
    print('sigma_z =\n', sigma_z)
    
    Sx1 = np.kron(sigma_x, np.diag([1., 0.]))
    Sx2 = np.kron(sigma_x, np.diag([0., 1.]))

    Sy1 = np.kron(sigma_y, np.diag([1., 0.]))
    Sy2 = np.kron(sigma_y, np.diag([0., 1.]))

    Sz1 = np.kron(sigma_z, np.diag([1., 0.]))
    Sz2 = np.kron(sigma_z, np.diag([0., 1.]))

    Sz = Sz1 + Sz2
    Sz_AF = Sz1 - Sz2
    
    print(Sz1)
    print(Sz2)
    print(Sz)

    # -- Sweep in interaction U
    
    #U_vec = np.arange(9., -0.5, -1.5)
    dU = -0.5
    U_vec = np.arange(9., 3.+dU, dU)
    
    res_AF = []
    res_FM = []

    M0 = np.zeros((4, 4))
    M_FM = -5.0 * Sz
    M_AF = -5.0 * Sz_AF
    mu_fm = 0.
    mu_af = 0.
    
    for U in U_vec:

        print('-'*72)
        print('U =', U)
        print('-'*72)

        mu0 = 0.5 * U # half-filling
        
        H_int = U * docc

        hs = HartreeFockSolver(
            e_k, beta, H_int=H_int, gf_struct=gf_struct)

        # -- Solve without seeding FM symmetry breaking
        print('-'*72)
        print('--> Restricted HF')
        hs.solve_newton(N_target=N_tot, M0=M0, mu0=mu0)

        print('-'*72)
        print('--> RPA response on non-magnetic HF solution')

        hr = HartreeFockResponse(hs)
        chi_FM = hr.response(Sz, Sz)
        chi_AF = hr.response(Sz_AF, Sz_AF)

        chi_m = hr.mode_decomposition()
                
        # -- Solve by seeding FM symmetry breaking

        print('-'*72)
        print('--> FM Unrestricted HF')
        print('M_FM =\n', M_FM)
        
        hs_fm = HartreeFockSolver(
            e_k, beta, H_int=H_int, gf_struct=gf_struct)
        
        hs_fm.solve_newton(N_target=N_tot, M0=M_FM, mu0=mu_fm)
        hs_fm.m_FM = hs_fm.expectation_value(Sz)
        hs_fm.m_AF = hs_fm.expectation_value(Sz_AF)
        hs_fm.U = U
        hs_fm.chi = chi_FM
        hs_fm.chi_m = chi_m

        mu_fm = hs_fm.mu
        M_FM = hs_fm.M.copy()

        print('-'*72)
        print('--> AF Unrestricted HF')
        print('M_AF =\n', M_AF)

        hs_af = HartreeFockSolver(
            e_k, beta, H_int=H_int, gf_struct=gf_struct)

        hs_af.solve_newton(N_target=N_tot, M0=M_AF, mu0=mu_af)
        hs_af.m_FM = hs_af.expectation_value(Sz)
        hs_af.m_AF = hs_af.expectation_value(Sz_AF)
        hs_af.U = U
        hs_af.chi = chi_AF
        
        mu_af = hs_af.mu
        M_AF = hs_af.M.copy()

        res_FM.append(hs_fm) 
        res_AF.append(hs_af) 

    # -- Check agreement between self-cons HF in FM and AF phase
    # -- and susceptibility of the non-spin polarized phase
    
    get_vec = lambda vec, attr : np.array([ getattr(el, attr) for el in vec])

    U = get_vec(res_AF, 'U')
    
    m_AF = get_vec(res_AF, 'm_AF')
    chi_AF = get_vec(res_AF, 'chi')
    
    assert( np.array_equal(np.abs(m_AF) < 1e-6, chi_AF.real > 0.) ), \
        "Error: The AF susceptibility and order parameter do not agree."

    m_FM = get_vec(res_FM, 'm_FM')
    chi_FM = get_vec(res_FM, 'chi')
    
    assert( np.array_equal(np.abs(m_FM) < 1e-6, chi_FM.real > 0.) ), \
        "Error: The AF susceptibility and order parameter do not agree."

    spl_chi_AF = InterpolatedUnivariateSpline(U[::-1], 1./chi_AF[::-1].real)
    U_AF = spl_chi_AF.roots()[0]
    print('U_AF =', U_AF)

    spl_chi_FM = InterpolatedUnivariateSpline(U[::-1], 1./chi_FM[::-1].real)
    U_FM = spl_chi_FM.roots()[0]
    print('U_FM =', U_FM)

    # -- Compare divergencies of chi with the
    
    chi_m = get_vec(res_FM, 'chi_m')
    
    e_AF = chi_m[:, -1] # AF is triply degenerate [spin-SU(2)]
    e_FM = chi_m[:, -4] # FM is the 4th eigen mode

    spl_e_AF = InterpolatedUnivariateSpline(U[::-1], 1./e_AF[::-1])
    U_AF_ref = spl_e_AF.roots()[0]
    print('U_AF_ref =', U_AF_ref)
    
    spl_e_FM = InterpolatedUnivariateSpline(U[::-1], 1./e_FM[::-1])
    U_FM_ref = spl_e_FM.roots()[0]
    print('U_FM_ref =', U_FM_ref)

    np.testing.assert_almost_equal(U_AF, U_AF_ref, decimal=4)
    np.testing.assert_almost_equal(U_FM, U_FM_ref, decimal=4)
    
    exit()

    import matplotlib.pyplot as plt

    plt.figure(figsize=(7, 7))

    for label, res in [('AF', res_AF), ('FM', res_FM)]:
        
        U = get_vec(res, 'U')
        m_AF = get_vec(res, 'm_AF').real
        m_FM = get_vec(res, 'm_FM').real
        chi = get_vec(res, 'chi').real

        E_tot = get_vec(res, 'E_tot').real
        E_kin = get_vec(res, 'E_kin').real
        E_int = get_vec(res, 'E_int').real

        Omega0 = get_vec(res, 'Omega0').real
        Omega = get_vec(res, 'Omega').real

        subp = [3, 2, 1]
        plt.subplot(*subp); subp[-1] += 1

        plt.title(
            r'$n=%2.2f$, $\beta=%2.2f$, $n_k=%s$' \
            % (N_tot, beta, str(n_k)))
        
        plt.plot(U, m_AF, '.-', label='%s $m_{AF}$' % label)
        plt.plot(U, m_FM, '.-', label='%s $m_{FM}$' % label)

        if label == 'AF':
            plt.plot([U_AF]*2, [0, np.max([np.max(m_AF), np.max(m_FM)])],
                     '.-r', label=r'$U_{c}^{AF}$')
            plt.plot([U_FM]*2, [0, np.max([np.max(m_AF), np.max(m_FM)])],
                     '.-r', label=r'$U_{c}^{FM}$')
        
        plt.xlabel(r'$U$')
        plt.ylabel(r'$m$')
        plt.legend(loc='best')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(U, 1./chi, '.-', label=r'$[\chi_{%s}]^{-1}$' % label)

        if label == 'AF':
            min_max = [(1./chi).min(), (1./chi).max()]
            plt.plot([U_AF]*2, min_max,
                     '.-r', label=r'$U_{c}^{AF}$')
            plt.plot([U_FM]*2, min_max,
                     '.-r', label=r'$U_{c}^{FM}$')

        plt.xlabel(r'$U$')
        plt.ylabel(r'$1/\chi$')
        plt.legend(loc='best')
        plt.grid(True)
            
        plt.subplot(*subp); subp[-1] += 1
        plt.plot(U, 0.5*E_tot, label=label + ' E_tot')
        plt.plot(U, 0.5*E_kin, label=label + ' E_kin')
        plt.plot(U, 0.5*E_int, label=label + ' E_int')
        plt.legend(loc='best', ncol=2)
        plt.xlabel(r'$U$')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(U, 0.5*Omega0, label=label + r' $\Omega_0$')
        plt.plot(U, 0.5*Omega, label=label + r' $\Omega$')
        plt.legend(loc='best', ncol=2)
        plt.xlabel(r'$U$')

        plt.subplot(*subp); subp[-1] += 1
        if hasattr(res[0], 'chi_m'):
            chi_m = get_vec(res, 'chi_m')
            plt.plot(U, 1./chi_m, '-', alpha=0.75)

        if label == 'AF':
            min_max = [(1./chi_m).min(), (1./chi_m).max()]
            plt.plot([U_AF]*2, min_max,
                     '.-r', label=r'$U_{c}^{AF}$')
            plt.plot([U_FM]*2, min_max,
                     '.-r', label=r'$U_{c}^{FM}$')
            
        plt.xlabel(r'$U$')
        plt.ylabel(r'Eigen values of $1-\chi^{(0)}_{AB} U_{AB}$')
        plt.grid(True)

    plt.tight_layout()
    plt.savefig('figure_1d_hubbard_hf_rpa.pdf')
    plt.show()

