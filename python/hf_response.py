
""" TRIQS, TPRF: Hartree-Fock static response function calculator

Author: Hugo U. R. Strand, hugo.strand@gmail.com (2018)
"""

# ----------------------------------------------------------------------

import itertools
import numpy as np

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
class HartreeResponse(BaseResponse):

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
        chi_ab = np.zeros((norb, norb), dtype=np.complex)
        for i1, i2 in itertools.product(range(norb), repeat=2):
            chi_ab[i1, i2] = chi_abcd[i1, i1, i2, i2]

        return chi_ab    

# ----------------------------------------------------------------------
