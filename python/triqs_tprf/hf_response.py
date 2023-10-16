# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2018 by The Simons Foundation
# Author: H. U.R. Strand
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

""" TRIQS, TPRF: Hartree-Fock static response function calculator

Author: Hugo U. R. Strand, hugo.strand@gmail.com (2018)
"""

# ----------------------------------------------------------------------

import sys
import itertools
import numpy as np

# ----------------------------------------------------------------------
class BaseResponse(object):

    def __init__(self, solver, eps=1e-9):

        print(self.logo())
        
        self.solver = solver
        self.eps = eps

        self.beta = self.solver.beta
        
        self.e_k = self.solver.e_k_MF.copy()
        self.n_k = len(self.e_k.mesh)
        
        self.norb = self.e_k.target_shape[0]
        self.shape_ab = self.e_k.target_shape
        self.shape_abcd = list(self.shape_ab) + list(self.shape_ab)
        self.shape_AB = [self.norb**2, self.norb**2]
        self.shape_kabcd = [self.n_k] + self.shape_abcd
        self.shape_kAB = [self.n_k] + self.shape_AB

        I_ab = np.eye(self.norb)
        self.e_k.data[:] -= self.solver.mu * I_ab

        print('norb =', self.norb)
        print('shape_abcd =', self.shape_abcd)
        print('shape_AB =', self.shape_AB)
        print('beta =', self.beta)
        
    # ------------------------------------------------------------------
    def _compute_drho_dop(self, op):

        drho_k = self._compute_drho_k_dop(op)
        drho = np.sum(drho_k, axis=0) / self.n_k

        return drho

    # ------------------------------------------------------------------
    def _compute_drho_k_dop(self, op):

        beta, eps = self.beta, self.eps
        
        E_p = self.e_k.data.copy() + eps*op[None, ...]
        E_m = self.e_k.data.copy() - eps*op[None, ...]

        e_p, U_p = np.linalg.eigh(E_p)
        e_m, U_m = np.linalg.eigh(E_m)

        fermi = lambda e, beta: 1./(np.exp(beta * e) + 1)

        rho_k_p = np.einsum('kab,kb,kcb->kac', U_p, fermi(e_p, beta), np.conj(U_p))
        rho_k_m = np.einsum('kab,kb,kcb->kac', U_m, fermi(e_m, beta), np.conj(U_m))

        drho_k = (rho_k_p - rho_k_m) / (2. * eps)

        return drho_k

    # ----------------------------------------------------------------------
    def _compute_chi0_ab(self):

        chi0_ab = np.zeros(self.shape_ab, dtype=complex)

        for a in range(self.norb):
            F_a = np.zeros(self.shape_ab)
            F_a[a, a] = 1.

            chi0_b = -np.diag(self._compute_drho_dop(F_a))
            chi0_ab[a] = chi0_b

        return chi0_ab

    # ----------------------------------------------------------------------
    def _compute_R_abcd(self, field_prefactor=1.):

        R_abcd = np.zeros(self.shape_abcd, dtype=complex)

        for a, b in itertools.product(list(range(self.norb)), repeat=2):

            F_ab = np.zeros(self.shape_ab, dtype=complex)
            F_ab[a, b] += field_prefactor
            F_ab[b, a] += np.conj(field_prefactor)

            R_cd = -self._compute_drho_dop(F_ab)
            R_abcd[b, a, :, :] = R_cd

        return R_abcd

    # ----------------------------------------------------------------------
    def _compute_R_kabcd(self, field_prefactor=1.):

        R_kabcd = np.zeros(self.shape_kabcd, dtype=complex)

        for a, b in itertools.product(list(range(self.norb)), repeat=2):

            F_ab = np.zeros(self.shape_ab, dtype=complex)
            F_ab[a, b] += field_prefactor
            F_ab[b, a] += np.conj(field_prefactor)

            R_kcd = -self._compute_drho_k_dop(F_ab)
            R_kabcd[:, b, a, :, :] = R_kcd

        return R_kabcd
    
    # ----------------------------------------------------------------------
    def _compute_chi0_abcd(self):

        chi0_kabcd = self._compute_chi0_kabcd()
        chi0_abcd = np.sum(chi0_kabcd, axis=0) / self.n_k

        return chi0_abcd

    # ----------------------------------------------------------------------
    def _compute_chi0_kabcd(self):

        R_r_kabcd = self._compute_R_kabcd(field_prefactor=1.0) 
        R_i_kabcd = self._compute_R_kabcd(field_prefactor=1.j)

        chi0_kabcd = 0.5 * (R_r_kabcd + R_i_kabcd.imag)

        return chi0_kabcd
    
    # ----------------------------------------------------------------------
    def logo(self):
        if 'UTF' in sys.stdout.encoding:
            txt = """
╔╦╗╦═╗╦╔═╗ ╔═╗  ┬ ┬┌─┐  ┬─┐┌─┐┌─┐
 ║ ╠╦╝║║═╬╗╚═╗  ├─┤├┤───├┬┘├─┘├─┤
 ╩ ╩╚═╩╚═╝╚╚═╝  ┴ ┴└    ┴└─┴  ┴ ┴
Triqs: Hartree-Fock Random Phase Approximation susceptibility
"""
        else:
            txt = r"""
 _____ ___ ___ ___  ___   _  _ ___    ___ ___  _
|_   _| _ \_ _/ _ \/ __| | || | __|__| _ \ _ \/_\
  | | |   /| | (_) \__ \ | __ | _|___|   /  _/ _ \
  |_| |_|_\___\__\_\___/ |_||_|_|    |_|_\_|/_/ \_\

Triqs: Hartree-Fock Random Phase Approximation susceptibility
"""
        return txt

# ----------------------------------------------------------------------
class HartreeFockResponse(BaseResponse):

    """ Hartree-Fock linear response calculator.

    Parameters
    ----------

    hartree_fock_solver : HartreeFockSolver instance
        Converged Hartree-Fock solver.

    eps : float
        Step size in finite difference linear response calculation.

    """
    
    def __init__(self, hartree_fock_solver, eps=1e-9):

        super(HartreeFockResponse, self).__init__(hartree_fock_solver)
        self.hfs = self.solver
        
        I_AB = np.matrix(np.eye(self.shape_AB[0]))
        U_AB = self._to_matrix_AB(self.hfs.U_abcd)
        chi0_AB = self._to_matrix_AB(self._compute_chi0_abcd())

        chi_AB = np.linalg.inv(I_AB - chi0_AB * U_AB) * chi0_AB

        self.chi0_abcd = self._to_tensor_abcd(chi0_AB)
        self.chi_abcd = self._to_tensor_abcd(chi_AB)

        #self.chi0_kabcd = self._compute_chi0_kabcd()
        #np.testing.assert_array_almost_equal(
        #    self.chi0_abcd, np.sum(self.chi0_kabcd, axis=0)/self.n_k)

    def mode_decomposition(self):
        
        U_AB = self._to_matrix_AB(self.hfs.U_abcd)
        chi0_AB = self._to_matrix_AB(self.chi0_abcd)

        e = np.linalg.eigvals(np.mat(chi0_AB) * np.mat(U_AB))
        
        idx = np.argsort(e.real)
        e = e[idx]
        e = 1./(1. - e.real)

        return e
        
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

    def compute_chi0_k(self):

        from triqs.gf import Gf

        chi0_k = Gf(mesh=self.e_k.mesh, target_shape=self.shape_abcd)
        chi0_k.data[:] = self.chi0_kabcd

        return chi0_k

# ----------------------------------------------------------------------
class HartreeResponse(BaseResponse):

    """ Hartree linear response calculator.

    Parameters
    ----------

    hartree_solver : HartreeSolver instance
        Converged Hartree solver.

    eps : float
        Step size in finite difference linear response calculation.

    """
    
    def __init__(self, hartree_solver, eps=1e-9):

        super(HartreeResponse, self).__init__(hartree_solver)
        
        I_ab = np.eye(self.norb)
        U_ab = np.mat(self.extract_dens_dens(self.solver.U_abcd))
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

    def extract_dens_dens(self, chi_abcd):
        norb = chi_abcd.shape[0]
        chi_ab = np.zeros((norb, norb), dtype=complex)
        for i1, i2 in itertools.product(list(range(norb)), repeat=2):
            chi_ab[i1, i2] = chi_abcd[i1, i1, i2, i2]

        return chi_ab    

# ----------------------------------------------------------------------
