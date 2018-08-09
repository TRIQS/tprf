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
def fermi_distribution_function(E, beta):
    I = np.eye(E.shape[0], dtype=np.complex)
    return np.mat(np.linalg.inv(expm(beta * E) + I))

# ----------------------------------------------------------------------
def local_density_matrix_g_wk(g_wk):

    # -- Local Gf
    g0_w = g0_wk[:, Idx(0,0,0)].copy()
    g0_w << 0.

    kmesh = g_wk.mesh.components[1]

    for k in kmesh:
        g0_w += g0_wk[:, k]
    g0_w /= len(kmesh)

    rho = g0_w.density()

    return rho

# ----------------------------------------------------------------------
def chi0_op_accurate(e_k, beta, op, eps=1e-9):

    shape = e_k.target_shape
    drho = np.zeros(shape, dtype=np.complex)

    for k in kmesh:

        e = e_k[k]

        rho_p = fermi_distribution_function(e + eps*op, beta)
        rho_m = fermi_distribution_function(e - eps*op, beta)
        drho += (rho_p - rho_m) / (2. * eps)
        
    drho /= len(kmesh)

    return -drho

# ----------------------------------------------------------------------
def rho_mf(e_k, beta, M=None, mu=None):

    if M is None:
        M = np.zeros(e_k.target_shape)
    if mu is None:
        mu = 0.

    I = np.eye(e_k.target_shape[0])
        
    E = e_k.data.copy() + M[None, ...] - mu * I[None, ...]

    e, V = np.linalg.eigh(E)

    fermi = lambda e, beta: 1./(np.exp(beta * e) + 1)

    rho = np.einsum('kab,kb,kcb->ac', V, fermi(e, beta), np.conj(V))
    rho /= len(e_k.mesh)

    return rho

# ----------------------------------------------------------------------
def fix_chemical_potential(e_k, beta, N, mu0=0.):

    def target_function(mu):
        norb = e_k.target_shape[0]
        rho = rho_mf(e_k, beta, -mu * np.eye(norb))
        n = np.sum(np.diag(rho)).real
        return n - N

    #from scipy.optimize import fsolve
    #mu = fsolve(target_function, mu0)

    from scipy.optimize import brentq
    mu = brentq(target_function, -10., 10.)
    #print 'fzero, result=', mu

    return [mu]

# ----------------------------------------------------------------------
def local_density_matrix(e_k, beta):

    shape = e_k.target_shape
    rho = np.zeros(shape, dtype=np.complex)

    for k in kmesh:
        rho += fermi_distribution_function(e_k[k], beta)
    rho /= len(kmesh)

    return rho

# ----------------------------------------------------------------------
def tr_op_chi_wk_op(op1, chi_wk, op2):
    
    chi_O1O2 = chi_wk[0, 0, 0, 0].copy()
    chi_O1O2.data[:] = np.einsum(
        'wqabcd,ab,cd->wq', chi_wk.data, op1, op2)[:, :]

    return chi_O1O2

# ----------------------------------------------------------------------
def U_matrix_kanamori_dens_dens(norb, U, J):

    Sz = np.kron(np.diag([+0.5, -0.5]), np.eye(norb/2))

    Oidx = np.diag(np.kron(np.eye(2), np.diag(range(norb/2))))
    Sidx = np.diag(Sz)

    U_mat = np.zeros((norb, norb))
    for i1, i2 in itertools.product(range(norb), repeat=2):
        s1, s2 = Sidx[i1], Sidx[i2]
        o1, o2 = Oidx[i1], Oidx[i2]

        if o1 == o2:
            if s1 != s2:
                U_mat[i1, i2] = U
        else:
            if s1 == s2:
                U_mat[i1, i2] = U - 3*J
            else:
                U_mat[i2, i1] = U - 2*J

    return U_mat

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
╔╦╗╦═╗╦╔═╗ ╔═╗  ┌┬┐┌─┐
 ║ ╠╦╝║║═╬╗╚═╗  │││├┤ 
 ╩ ╩╚═╩╚═╝╚╚═╝  ┴ ┴└  
TRIQS: Mean-Field solver
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
            assert( H_int is None )

        if H_int is not None:
            assert( gf_struct is not None )

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
    def update_density_matrix(self):

        e, V = np.linalg.eigh(self.e_k_MF.data)
        e -= self.mu
        
        fermi = lambda e : 1./(np.exp(self.beta * e) + 1)

        rho = np.einsum('kab,kb,kcb->ac', V, fermi(e), np.conj(V))
        rho /= len(self.e_k_MF.mesh)

        self.rho = rho
        self.N_tot = np.sum(np.diag(self.rho))
        
        return np.diag(self.rho).real

    # ------------------------------------------------------------------
    def density_matrix_step(self, rho_diag):
        
        self.update_mean_field(rho_diag)
        self.update_mean_field_dispersion()
        self.update_chemical_potential(self.N_target, mu0=self.mu)
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
            rho_diag_new = self.density_matrix_step(rho_diag)

            drho = np.linalg.norm(rho_diag_old - rho_diag_new) / np.linalg.norm(rho_diag_old)
            print 'MF: iter, drho = %3i, %2.2E' % (idx, drho)
            
            if drho < tol:
                print 'MF: Converged drho = %3.3E\n' % drho
                print self.__str__()
                break

            rho_diag = (1. - mixing) * rho_diag_old + mixing * rho_diag_new

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
            rho_diag_new = self.density_matrix_step(rho_diag)
            return rho_diag_new - rho_diag
            
        rho_diag = fsolve(target_function, rho0_diag)

        print self.__str__()
        
        self.density_matrix_step(rho_diag)

    # ------------------------------------------------------------------
    def __str__(self):

        txt = 'MF: Solver state\n' + \
            'beta = ' + str(self.beta) + '\n' + \
            'N_target = ' + str(self.N_tot) + '\n' + \
            'N_tot    = ' + str(np.sum(np.diag(self.rho.real))) + '\n' + \
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
class HartreeResponse(BaseResponse):

    def __init__(self, hartree_solver, eps=1e-9):

        super(HartreeResponse, self).__init__(hartree_solver)
        
        I_ab = np.eye(self.norb)
        U_ab = np.mat(self.solver.U_ab)
        chi0_ab = np.mat(self._compute_chi0_ab())
        chi_ab = chi0_ab * np.linalg.inv(I_ab - U_ab * chi0_ab)

        self.chi0_ab = np.array(chi0_ab)
        self.chi_ab = np.array(chi_ab)

    def __check_op(self, op):
        """ Operators have to be diagonal in the Hartree approx """
        
        assert( op.shape == self.e_k.target_shape )
        np.testing.assert_almost_equal(
            op - np.diag(np.diag(op)), np.zeros_like(op))
        
    def bare_response(self, op1, op2):

        self.__check_op(op1)
        self.__check_op(op2)
        
        chi0_op1op2 = np.einsum('aa,ab,bb->', op1, self.chi0_ab, op2)

        return chi0_op1op2

    def response(self, op1, op2):

        self.__check_op(op1)
        self.__check_op(op2)
        
        chi_op1op2 = np.einsum('aa,ab,bb->', op1, self.chi_ab, op2)

        return chi_op1op2        

# ----------------------------------------------------------------------
class HartreeFockResponse(BaseResponse):

    def __init__(self, hartree_fock_solver, eps=1e-9):

        super(HartreeFockResponse, self).__init__(hartree_fock_solver)
        self.hfs = self.solver
        
        I_AB = np.matrix(np.eye(self.shape_AB[0]))
        U_AB = self._to_matrix_AB(self.hfs.U_abcd)
        chi0_AB = self._to_matrix_AB(self._compute_chi0_abcd())
        
        chi_AB = np.linalg.inv(I_AB - chi0_AB * U_AB) * chi0_AB

        self.chi0_abcd = self._to_tensor_abcd(chi0_AB)
        self.chi_abcd = self._to_tensor_abcd(chi_AB)

    def _to_matrix_AB(self, tensor_abcd):
        matrix_AB = np.matrix(tensor_abcd.reshape(self.shape_AB))
        return matrix_AB

    def _to_tensor_abcd(self, matrix_AB):
        tensor_abcd = np.array(matrix_AB).reshape(self.shape_abcd)
        return tensor_abcd
    
    def __check_op(self, op):
        assert( op.shape == self.shape_ab )
        
    def bare_response(self, op1, op2):

        self.__check_op(op1)
        self.__check_op(op2)
        
        chi0_op1op2 = np.einsum('ab,abcd,cd->', op1, self.chi0_abcd, op2)

        return chi0_op1op2

    def response(self, op1, op2):

        self.__check_op(op1)
        self.__check_op(op2)
        
        chi_op1op2 = np.einsum('ab,abcd,cd->', op1, self.chi_abcd, op2)

        return chi_op1op2

# ----------------------------------------------------------------------
