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
        self.n_k = len(self.e_k.mesh)

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
    def update_mean_field(self, rho_ab):

        self.M = np.einsum('abcd,cd->ba', -self.U_abcd, rho_ab)
        return self.M
    
    # ------------------------------------------------------------------
    def update_mean_field_dispersion(self):
        self.e_k_MF.data[:] = self.e_k.data + self.M[None, ...]

    # ------------------------------------------------------------------
    def update_chemical_potential(self, N_target, mu0=None):

        if mu0 is None:
            mu0 = self.mu
            
        e = np.linalg.eigvalsh(self.e_k_MF.data)

        fermi = lambda e : 1./(np.exp(self.beta * e) + 1)

        def target_function(mu):
            n = np.sum(fermi(e - mu)) / self.n_k
            return n - N_target

        mu = brentq(target_function, self.mu_min, self.mu_max)

        self.mu = mu        

    # ------------------------------------------------------------------
    def update_momentum_density_matrix(self):

        e, V = np.linalg.eigh(self.e_k_MF.data)
        e -= self.mu
        
        fermi = lambda e : 1./(np.exp(self.beta * e) + 1)
        self.rho_kab = np.einsum('kab,kb,kcb->kac', V, fermi(e), np.conj(V))

        return self.rho_kab

    # ------------------------------------------------------------------
    def update_density_matrix(self):

        rho_kab = self.update_momentum_density_matrix()
        self.rho_ab = np.sum(rho_kab, axis=0) / self.n_k
        self.N_tot = np.sum(np.diag(self.rho_ab))
        
        return self.rho_ab

    # ------------------------------------------------------------------
    def update_non_int_free_energy(self):

        e = np.linalg.eigvalsh(self.e_k_MF.data)
        e -= self.mu
        
        self.Omega0 = -1./self.beta * np.sum( np.log(1. + np.exp(-self.beta*e)) )
        self.Omega0 /= len(self.e_k_MF.mesh)
        
        return self.Omega0

    # ------------------------------------------------------------------
    def update_kinetic_energy(self):

        rho_kab = self.update_momentum_density_matrix()
        self.E_kin = np.einsum('kab,kba->', self.e_k.data, rho_kab)
        self.E_kin /= len(self.e_k.mesh)

        return self.E_kin
        
    # ------------------------------------------------------------------
    def update_interaction_energy(self):

        #self.E_int = 0.5 * np.einsum('aa,ab,bb->', self.rho, self.U_ab, self.rho)
        self.E_int = 0.5 * np.einsum(
            'ab,abcd,cd->', self.rho_ab, self.U_abcd, self.rho_ab)
        return self.E_int        
        
    # ------------------------------------------------------------------
    def update_total_energy(self):

        E_kin = self.update_kinetic_energy()
        E_int = self.update_interaction_energy()
        self.E_tot = E_kin + E_int

        Omega0 = self.update_non_int_free_energy()
        self.Omega = Omega0 + E_int

        return self.E_tot
        
    # ------------------------------------------------------------------
    def density_matrix_step(self, rho_vec, N_target=None):

        rho_ab = self.vec2mat(rho_vec)
        
        self.update_mean_field(rho_ab)
        self.update_mean_field_dispersion()

        if N_target is not None:
            self.update_chemical_potential(N_target, mu0=self.mu)
            
        rho_ab = self.update_density_matrix()

        rho_vec = self.mat2vec(rho_ab)

        return rho_vec

    # ------------------------------------------------------------------
    def solve_setup(self, N_target, M0=None, mu0=None):
    
        self.N_target = N_target

        if mu0 is None:
            self.update_chemical_potential(self.N_target, mu0=0.)
        else:
            self.mu = mu0

        if M0 is None:
            M0 = np.zeros(self.target_shape)

        M0 = np.array(M0, dtype=np.complex)

        self.M = np.copy(M0)
        self.update_mean_field_dispersion()
        self.update_chemical_potential(self.N_target)
        rho0_ab = self.update_density_matrix()

        rho0_vec = self.mat2vec(rho0_ab)
        
        return rho0_vec

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
        
        rho_vec = self.solve_setup(N_target, M0, mu0)

        rho_iter = []
        
        for idx in xrange(self.nitermax):

            rho_iter.append(rho_vec)

            rho_vec_old = np.copy(rho_vec)
            rho_vec_new = self.density_matrix_step(rho_vec, N_target)

            norm = np.linalg.norm(rho_vec_old)
            drho = np.linalg.norm(rho_vec_old - rho_vec_new) / norm
            
            print 'MF: iter, drho = %3i, %2.2E' % (idx, drho)
            
            if drho < tol:
                print 'MF: Converged drho = %3.3E\n' % drho
                break

            rho_vec = (1. - mixing) * rho_vec_old + mixing * rho_vec_new

        self.update_total_energy()
        print self.__str__()

        rho_iter = np.array(rho_iter)
        return rho_iter

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
        
        rho0_vec = self.solve_setup(N_target, M0, mu0)
        
        def target_function(rho_vec):
            rho_vec_new = self.density_matrix_step(rho_vec, N_target)
            return rho_vec_new - rho_vec
            
        rho_vec = fsolve(target_function, rho0_vec)

        self.update_total_energy()
        print self.__str__()
        
        self.density_matrix_step(rho_vec)

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

        if M0 is None:
            M0 = np.zeros(self.target_shape)

        M0 = np.array(M0, dtype=np.complex)
            
        self.M = np.copy(M0)
        self.update_mean_field_dispersion()
        rho0_ab = self.update_density_matrix()
        rho0_vec = self.mat2vec(rho0_ab)
        
        def target_function(rho_vec):
            rho_vec_new = self.density_matrix_step(rho_vec)
            return rho_vec_new - rho_vec
            
        rho_vec = fsolve(target_function, rho0_vec)

        self.update_total_energy()
        print self.__str__()
        
        self.density_matrix_step(rho_vec)
        
    # ------------------------------------------------------------------
    def total_density(self):
        return self.N_tot

    # ------------------------------------------------------------------
    def density_matrix(self):
        return self.rho_ab

    # ------------------------------------------------------------------
    def mean_field_matrix(self):
        return self.M

    # ------------------------------------------------------------------
    def chemical_potential(self):
        return self.mu

    # ------------------------------------------------------------------
    def expectation_value(self, op_ab):
        return np.einsum('ab,ba->', op_ab, self.rho_ab)
    
    # ------------------------------------------------------------------
    def __str__(self):

        txt = 'MF: Solver state\n' + \
            'beta = ' + str(self.beta) + '\n' + \
            'N_tot = ' + str(self.N_tot) + '\n' + \
            'E_tot    = ' + str(self.E_tot) + '\n' + \
            'E_int    = ' + str(self.E_int) + '\n' + \
            'E_kin    = ' + str(self.E_kin) + '\n' + \
            'Omega0   = ' + str(self.Omega0) + '\n' + \
            'rho_ab =\n' + str(self.rho_ab) + '\n' + \
            'mu = ' + str(self.mu) + '\n' + \
            'M =\n' + str(self.M) + '\n' + \
            'M - mu =\n' + str(self.M - self.mu * np.eye(self.norb)) + '\n'

        return txt

    # ------------------------------------------------------------------
    def mat2vec(self, mat):
        assert( len(mat.shape) == 2 )
        return np.diag(mat).real

    # ------------------------------------------------------------------
    def vec2mat(self, vec):
        assert( len(vec.shape) == 1 )
        return np.diag(vec.real)    
    
# ----------------------------------------------------------------------
class HartreeFockSolver(HartreeSolver):

    # ------------------------------------------------------------------
    def mat2vec(self, mat):
        assert( len(mat.shape) == 2 )

        diag = np.diag(mat).real

        idx = np.triu_indices(mat.shape[0], k=1)
        up_tri = mat[idx]

        vec = np.concatenate((diag, up_tri.real, up_tri.imag))        

        if False:
            print '-'*72
            print 'vec =', vec
            print 'mat =\n', mat
            print '-'*72
        
        return vec

    # ------------------------------------------------------------------
    def vec2mat(self, vec):

        assert( len(vec.shape) == 1 )

        n = self.norb
        diag, up_tri = vec[:n], vec[n:]

        mat = np.zeros((n, n), dtype=np.complex)

        re, im = np.split(up_tri, 2)

        idx = np.triu_indices(mat.shape[0], k=1)
        mat[idx] = re + 1.j * im

        mat = mat + np.conj(mat).T
        mat += np.diag(diag)

        if False:
            print '-'*72
            print 'vec =', vec
            print 'mat =\n', mat
            print '-'*72
        
        return mat
    
# ----------------------------------------------------------------------
