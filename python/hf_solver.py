# -*- coding: utf-8 -*-

""" TRIQS: Density-Density Mean-Field solver

Author: Hugo U. R. Strand, hugo.strand@gmail.com (2018)
"""

# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from scipy.optimize import fsolve
from scipy.optimize import brentq

# ----------------------------------------------------------------------

from triqs_tprf.rpa_tensor import get_rpa_tensor
from triqs_tprf.rpa_tensor import fundamental_operators_from_gf_struct

# ----------------------------------------------------------------------
def extract_dens_dens(chi_abcd):
    norb = chi_abcd.shape[0]
    chi_ab = np.zeros((norb, norb), dtype=np.complex)
    for i1, i2 in itertools.product(range(norb), repeat=2):
        chi_ab[i1, i2] = chi_abcd[i1, i1, i2, i2]

    return chi_ab

# ----------------------------------------------------------------------
class HartreeSolver(object):

    # ------------------------------------------------------------------
    def __init__(self, e_k, beta, H_int=None, gf_struct=None,
                 mu0=0., mu_max=10, mu_min=-10.):
        """ TRIQS: Hartree solver (density-density)

        Parameters:

        e_k : single-particle dispersion
        U_mat : density-density interaction matrix
        mu_min, mu_max : range for chemical potential search
        
        """

        logo = """
╔╦╗╦═╗╦╔═╗ ╔═╗  ┬ ┬┌─┐
 ║ ╠╦╝║║═╬╗╚═╗  ├─┤├┤ 
 ╩ ╩╚═╩╚═╝╚╚═╝  ┴ ┴└  
TRIQS: Hartree-Fock solver
"""
        print logo
        
        self.mu = mu0
        self.beta = beta
        self.mu_max = mu_max
        self.mu_min = mu_min
        
        self.e_k = e_k.copy()
        self.e_k_MF = e_k.copy()

        self.target_shape = self.e_k.target_shape

        self.shape_ab = self.e_k.target_shape
        self.shape_abcd = list(self.shape_ab) * 2

        self.norb = self.target_shape[0]

        print 'bands =', self.norb
        print 'n_k =', len(self.e_k.mesh)
        print

        if gf_struct is None:
            assert( H_int is None ), \
                'Error: gf_struct = None, but H_int is not None'

        if H_int is not None:
            assert( gf_struct is not None ), \
                'Error: H_int = None, but gf_struct is not None'

            fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)
            self.U_abcd = get_rpa_tensor(H_int, fundamental_operators)
            self.U_ab = extract_dens_dens(self.U_abcd)

        else:
            self.U_abcd = np.zeros(self.shape_abcd)
            self.U_ab = np.zeros(self.shape_ab)

        assert( self.e_k.target_shape == self.U_ab.shape )
        
    # ------------------------------------------------------------------
    def update_mean_field(self, rho_diag):
        self.M = np.diag( np.dot( -self.U_ab, rho_diag) )
        return self.M
    
    # ------------------------------------------------------------------
    def update_mean_field_dispersion(self):
        self.e_k_MF.data[:] = self.e_k.data + self.M[None, ...]

    # ------------------------------------------------------------------
    def update_chemical_potential(self, N_target, mu0=None):

        if mu0 is None:
            mu0 = self.mu
            
        n_k = len(self.e_k_MF.mesh)
        e = np.linalg.eigvalsh(self.e_k_MF.data)

        fermi = lambda e : 1./(np.exp(self.beta * e) + 1)

        def target_function(mu):
            n = np.sum(fermi(e - mu)) / n_k
            return n - N_target

        mu = brentq(target_function, self.mu_min, self.mu_max)

        self.mu = mu        

    # ------------------------------------------------------------------
    def update_momentum_density_matrix(self):

        e, V = np.linalg.eigh(self.e_k_MF.data)
        e -= self.mu
        
        fermi = lambda e : 1./(np.exp(self.beta * e) + 1)

        self.rho_k = np.einsum('kab,kb,kcb->kac', V, fermi(e), np.conj(V))

        return self.rho_k

    # ------------------------------------------------------------------
    def update_density_matrix(self):

        rho_k = self.update_momentum_density_matrix()
        rho = np.sum(rho_k, axis=0) / len(self.e_k_MF.mesh)

        self.rho = rho
        self.N_tot = np.sum(np.diag(self.rho))
        
        return np.diag(self.rho).real

    # ------------------------------------------------------------------
    def update_non_int_free_energy(self):

        e = np.linalg.eigvalsh(self.e_k_MF.data)
        e -= self.mu
        
        self.Omega0 = -1./self.beta * np.sum( np.log(1. + np.exp(-self.beta*e)) )
        self.Omega0 /= len(self.e_k_MF.mesh)
        
        return self.Omega0

    # ------------------------------------------------------------------
    def update_kinetic_energy(self):

        rho_k = self.update_momentum_density_matrix()
        self.E_kin = np.einsum('kab,kba->', self.e_k.data, self.rho_k)
        self.E_kin /= len(self.e_k.mesh)

        return self.E_kin
        
    # ------------------------------------------------------------------
    def update_interaction_energy(self):

        self.E_int = 0.5 * np.einsum('aa,ab,bb->', self.rho, self.U_ab, self.rho)

        return self.E_int        
        
    # ------------------------------------------------------------------
    def update_total_energy(self):
        E_kin = self.update_kinetic_energy()
        E_int = self.update_interaction_energy()
        self.E_tot = E_kin + E_int

        self.update_non_int_free_energy()
        self.Omega = self.Omega0 + self.E_int

        return self.E_tot
        
    # ------------------------------------------------------------------
    def density_matrix_step(self, rho_diag, N_target=None):
        
        self.update_mean_field(rho_diag)
        self.update_mean_field_dispersion()

        if N_target is not None:
            self.update_chemical_potential(N_target, mu0=self.mu)
            
        rho_diag = self.update_density_matrix()

        return rho_diag

    # ------------------------------------------------------------------
    def solve_setup(self, N_target, M0=None, mu0=None):
    
        self.N_target = N_target

        if mu0 is None:
            self.update_chemical_potential(self.N_target, mu0=0.)
        else:
            self.mu = mu0

        if M0 is None:
            M0 = np.zeros(self.target_shape)

        self.M = np.copy(M0)
        self.update_mean_field_dispersion()
        self.update_chemical_potential(self.N_target)
        rho0_diag = self.update_density_matrix()

        return rho0_diag

    # ------------------------------------------------------------------
    def solve_iter(self, N_target, M0=None, mu0=None,
                   nitermax=100, mixing=0.5, tol=1e-9):
        """
        Parameters:

        beta : Inverse temperature
        N_tot : total density
        M0 : Initial mean-field (0 if None)
        mu0 : Initial chemical potential

        nitermax : maximal number of self consistent iterations
        mixing : linear mixing parameter
        tol : convergence in relative change of the density matrix

        """

        print 'MF: Forward iteration'
        print 'nitermax =', nitermax
        print 'mixing =', mixing
        print 'tol =', tol
        print
        
        assert( mixing >= 0. )
        assert( mixing <= 1. )

        self.mixing = mixing
        self.nitermax = nitermax
        
        rho_diag = self.solve_setup(N_target, M0, mu0)

        rho_vec = []
        
        for idx in xrange(self.nitermax):

            rho_vec.append(rho_diag)

            rho_diag_old = np.copy(rho_diag)
            rho_diag_new = self.density_matrix_step(rho_diag, N_target)

            drho = np.linalg.norm(rho_diag_old - rho_diag_new) / np.linalg.norm(rho_diag_old)
            print 'MF: iter, drho = %3i, %2.2E' % (idx, drho)
            
            if drho < tol:
                print 'MF: Converged drho = %3.3E\n' % drho
                print self.__str__()
                break

            rho_diag = (1. - mixing) * rho_diag_old + mixing * rho_diag_new

        self.update_total_energy()
        print self.__str__()

        rho_vec = np.array(rho_vec)
        return rho_vec

    # ------------------------------------------------------------------
    def solve_newton(self, N_target, M0=None, mu0=None):
        """
        Parameters:

        beta : Inverse temperature
        N_tot : total density
        M0 : Initial mean-field (0 if None)
        mu0 : Initial chemical potential
        """

        print 'MF: Newton solver'
        print
        
        rho0_diag = self.solve_setup(N_target, M0, mu0)
        
        def target_function(rho_diag):
            rho_diag_new = self.density_matrix_step(rho_diag, N_target)
            return rho_diag_new - rho_diag
            
        rho_diag = fsolve(target_function, rho0_diag)

        self.update_total_energy()
        print self.__str__()
        
        self.density_matrix_step(rho_diag)

    # ------------------------------------------------------------------
    def solve_newton_mu(self, mu, M0=None):
        """
        Parameters:

        beta : Inverse temperature
        N_tot : total density
        M0 : Initial mean-field (0 if None)
        mu0 : Initial chemical potential
        """

        print 'MF: Newton solver'
        print

        self.mu = mu
        #rho0_diag = self.solve_setup(N_target, M0, mu0)

        if M0 is None:
            M0 = np.zeros(self.target_shape)
            
        self.M = np.copy(M0)
        self.update_mean_field_dispersion()
        rho0_diag = self.update_density_matrix()
        
        def target_function(rho_diag):
            rho_diag_new = self.density_matrix_step(rho_diag)
            return rho_diag_new - rho_diag
            
        rho_diag = fsolve(target_function, rho0_diag)

        self.update_total_energy()
        print self.__str__()
        
        self.density_matrix_step(rho_diag)
        
    # ------------------------------------------------------------------
    def __str__(self):

        txt = 'MF: Solver state\n' + \
            'beta = ' + str(self.beta) + '\n' + \
            'N_target = ' + str(self.N_tot) + '\n' + \
            'N_tot    = ' + str(np.sum(np.diag(self.rho.real))) + '\n' + \
            'E_tot    = ' + str(self.E_tot) + '\n' + \
            'E_int    = ' + str(self.E_int) + '\n' + \
            'E_kin    = ' + str(self.E_kin) + '\n' + \
            'Omega0   = ' + str(self.Omega0) + '\n' + \
            'rho =\n' + str(self.rho) + '\n' + \
            'mu = ' + str(self.mu) + '\n' + \
            'M =\n' + str(self.M) + '\n' + \
            'M - mu =\n' + str(self.M - self.mu * np.eye(self.norb)) + '\n'

        return txt
        
# ----------------------------------------------------------------------
class BaseResponse(object):

    def __init__(self, solver, eps=1e-9):

        self.solver = solver
        self.eps = eps

        self.beta = self.solver.beta
        
        self.e_k = self.solver.e_k_MF.copy()
        
        self.norb = self.e_k.target_shape[0]
        self.shape_ab = self.e_k.target_shape
        self.shape_abcd = list(self.shape_ab) + list(self.shape_ab)
        self.shape_AB = (self.norb**2, self.norb**2)

        I_ab = np.eye(self.norb)
        self.e_k.data[:] -= self.solver.mu * I_ab
        
    # ------------------------------------------------------------------
    def _compute_drho_dop(self, op):

        beta, eps = self.beta, self.eps
        
        E_p = self.e_k.data.copy() + eps*op[None, ...]
        E_m = self.e_k.data.copy() - eps*op[None, ...]

        e_p, U_p = np.linalg.eigh(E_p)
        e_m, U_m = np.linalg.eigh(E_m)

        fermi = lambda e, beta: 1./(np.exp(beta * e) + 1)

        rho_p = np.einsum('kab,kb,kcb->kac', U_p, fermi(e_p, beta), np.conj(U_p))
        rho_m = np.einsum('kab,kb,kcb->kac', U_m, fermi(e_m, beta), np.conj(U_m))

        drho = (rho_p - rho_m) / (2. * eps)
        drho = np.sum(drho, axis=0) / len(self.e_k.mesh)

        return drho

    # ----------------------------------------------------------------------
    def _compute_chi0_ab(self):

        chi0_ab = np.zeros(self.shape_ab, dtype=np.complex)

        for a in range(self.norb):
            F_a = np.zeros(self.shape_ab)
            F_a[a, a] = 1.

            chi0_b = -np.diag(self._compute_drho_dop(F_a))
            chi0_ab[a] = chi0_b

        return chi0_ab

    # ----------------------------------------------------------------------
    def _compute_R_abcd(self, field_prefactor=1.):

        R_abcd = np.zeros(self.shape_abcd, dtype=np.complex)

        for a, b in itertools.product(range(self.norb), repeat=2):

            F_ab = np.zeros(self.shape_ab, dtype=np.complex)
            F_ab[a, b] += field_prefactor
            F_ab[b, a] += np.conj(field_prefactor)

            R_cd = -self._compute_drho_dop(F_ab)
            R_abcd[b, a, :, :] = R_cd

        return R_abcd

    # ----------------------------------------------------------------------
    def _compute_chi0_abcd(self):

        R_r_abcd = self._compute_R_abcd(field_prefactor=1.0) 
        R_i_abcd = self._compute_R_abcd(field_prefactor=1.j)

        chi0_abcd = 0.5 * (R_r_abcd + R_i_abcd.imag)

        return chi0_abcd

# ----------------------------------------------------------------------
