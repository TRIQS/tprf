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

""" TRIQS: Hartree-Fock Solver

for systems with general dispersion and local interaction

Author: Hugo U. R. Strand, hugo.strand@gmail.com (2018)
"""

# ----------------------------------------------------------------------

import sys
import itertools
import numpy as np

# ----------------------------------------------------------------------

from scipy.optimize import fsolve
from scipy.optimize import brentq

# ----------------------------------------------------------------------

from triqs.gf.block_gf import fix_gf_struct_type
import triqs.utility.mpi as mpi

# ----------------------------------------------------------------------

from triqs_tprf.rpa_tensor import get_rpa_tensor
from triqs_tprf.rpa_tensor import fundamental_operators_from_gf_struct
from triqs_tprf.OperatorUtils import is_operator_composed_of_only_fundamental_operators

# ----------------------------------------------------------------------
class HartreeFockSolver(object):

    """ Hartree-Fock solver for local interactions.

    Parameters
    ----------

    e_k : TRIQS Greens function on a Brillouin zone mesh
        Single-particle dispersion.

    beta : float
        inverse temperature

    H_int : TRIQS Operator instance
        Local interaction Hamiltonian

    gf_struct : list of pairs of block index and its size
        gf_struct fixing orbital order between e_k and H_int

    mu0 : float
        chemical potential

    mu_min, mu_max : float, optional
        range for chemical potential search.

    """
        
    # ------------------------------------------------------------------
    def __init__(self, e_k, beta, H_int=None, gf_struct=None,
                 mu0=0., mu_max=None, mu_min=None):

        if mpi.is_master_node():
            print(self.logo())

        gf_struct = fix_gf_struct_type(gf_struct)
        
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
        self.triu_idxs = np.triu_indices(self.norb, k=1)

        if mpi.is_master_node():
            print('beta =', self.beta)
            print(f'mu = {self.mu} -- (min/max = {self.mu_min}/{self.mu_max})')
            print('bands =', self.norb)
            print('n_k =', len(self.e_k.mesh))
            print('H_int =', H_int)
            print()

        if gf_struct is None:
            assert( H_int is None ), \
                'Error: gf_struct = None, but H_int is not None'

        if H_int is not None:
            assert( gf_struct is not None ), \
                'Error: H_int = None, but gf_struct is not None'

            fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)

            if mpi.is_master_node():
                print('gf_struct =', gf_struct)
                print('fundamental_operators =', fundamental_operators)
                print()

            assert( is_operator_composed_of_only_fundamental_operators(
                H_int, fundamental_operators) ), \
                'Error: H_int is incompatible with gf_struct and its fundamental_operators'
            
            self.U_abcd = get_rpa_tensor(H_int, fundamental_operators)

        else:
            self.U_abcd = np.zeros(self.shape_abcd)

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

        mu_min = self.mu_min if self.mu_min is not None else e.min()
        mu_max = self.mu_max if self.mu_max is not None else e.max()

        mu = brentq(target_function, mu_min, mu_max)

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

        M0 = np.array(M0, dtype=complex)

        self.M = np.copy(M0)
        self.update_mean_field_dispersion()
        self.update_chemical_potential(self.N_target)
        rho0_ab = self.update_density_matrix()

        rho0_vec = self.mat2vec(rho0_ab)
        
        return rho0_vec

    # ------------------------------------------------------------------
    def solve_iter(self, N_target, M0=None, mu0=None,
                   nitermax=100, mixing=0.5, tol=1e-9):
        """ Solve the HF-equations using forward recursion at fixed density.     

        Parameters
        ----------

        N_target : float
            Total density.

        M0 : ndarray (2D), optional
            Initial mean-field (0 if None).

        mu0 : float, optional
            Initial chemical potential.

        nitermax : int, optional
            Maximal number of self consistent iterations.

        mixing : float, optional
            Linear mixing parameter.

        tol : float, optional
            Convergence in relative change of the density matrix.

        Returns
        -------

        rho : ndarray (2D)
            Local density matrix.

        """

        print('MF: Forward iteration')
        print('nitermax =', nitermax)
        print('mixing =', mixing)
        print('tol =', tol)
        print()
        
        assert( mixing >= 0. )
        assert( mixing <= 1. )

        self.mixing = mixing
        self.nitermax = nitermax
        
        rho_vec = self.solve_setup(N_target, M0, mu0)

        rho_iter = []
        
        for idx in range(self.nitermax):

            rho_iter.append(rho_vec)

            rho_vec_old = np.copy(rho_vec)
            rho_vec_new = self.density_matrix_step(rho_vec, N_target)

            norm = np.linalg.norm(rho_vec_old)
            drho = np.linalg.norm(rho_vec_old - rho_vec_new) / norm
            
            print('MF: iter, drho = %3i, %2.2E' % (idx, drho))
            
            if drho < tol:
                print('MF: Converged drho = %3.3E\n' % drho)
                break

            rho_vec = (1. - mixing) * rho_vec_old + mixing * rho_vec_new

        self.update_total_energy()
        print(self.__str__())

        rho_iter = np.array(rho_iter)
        return rho_iter

    # ------------------------------------------------------------------
    def solve_newton(self, N_target, M0=None, mu0=None):

        """ Solve the HF-equations using a quasi Newton method at fixed density.

        Parameters
        ----------

        N_tot : float
            Total density.

        M0 : ndarray (2D), optional
            Initial mean-field (0 if None).

        mu0 : float, optional
            Initial chemical potential.

        """

        print('MF: Newton solver')
        print()
        
        rho0_vec = self.solve_setup(N_target, M0, mu0)
        
        def target_function(rho_vec):
            rho_vec_new = self.density_matrix_step(rho_vec, N_target)
            return rho_vec_new - rho_vec
            
        rho_vec = fsolve(target_function, rho0_vec)

        self.update_total_energy()
        print(self.__str__())
        
        self.density_matrix_step(rho_vec)

    # ------------------------------------------------------------------
    def solve_newton_mu(self, mu, M0=None):

        """ Solve the HF-equations using a quasi Newton method at fixed chemical potential.

        Parameters
        ----------

        mu0 : float
            Initial chemical potential.

        M0 : ndarray (2D), optional
            Initial mean-field (0 if None).

        """

        print('MF: Newton solver')
        print()

        self.mu = mu

        if M0 is None:
            M0 = np.zeros(self.target_shape)

        M0 = np.array(M0, dtype=complex)
            
        self.M = np.copy(M0)
        self.update_mean_field_dispersion()
        rho0_ab = self.update_density_matrix()
        rho0_vec = self.mat2vec(rho0_ab)
        
        def target_function(rho_vec):
            rho_vec_new = self.density_matrix_step(rho_vec)
            return rho_vec_new - rho_vec
            
        rho_vec = fsolve(target_function, rho0_vec)

        self.update_total_energy()
        print(self.__str__())
        
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
            'Omega    = ' + str(self.Omega) + '\n' + \
            'rho_ab =\n' + str(self.rho_ab) + '\n' + \
            'mu = ' + str(self.mu) + '\n' + \
            'M =\n' + str(self.M) + '\n' + \
            'M - mu =\n' + str(self.M - self.mu * np.eye(self.norb)) + '\n'

        return txt

    # ------------------------------------------------------------------
    def mat2vec(self, mat):
        
        r""" Converts a unitary matrix to a vector representation
        with the order

        1. the real diagonal entries
        2. the real part of the upper triangular entries
        3. the imaginary part of the upper triangular entries

        """
        
        assert( len(mat.shape) == 2 )

        diag = np.diag(mat).real
        up_tri = mat[self.triu_idxs]

        vec = np.concatenate((diag, up_tri.real, up_tri.imag))        

        return vec

    # ------------------------------------------------------------------
    def vec2mat(self, vec):
        
        """ Converts from vector representation to a unitary matrix
        see mat2vec(...) for details. """

        assert( len(vec.shape) == 1 )

        mat = np.zeros(self.shape_ab, dtype=complex)

        diag, up_tri = vec[:self.norb], vec[self.norb:]
        re, im = np.split(up_tri, 2)
        
        mat[self.triu_idxs] = re + 1.j * im

        mat += np.conj(mat).T # reconstruct lower triangle
        mat += np.diag(diag) # add diagonal

        return mat
    
    # ------------------------------------------------------------------
    def logo(self):
        if 'UTF' in sys.stdout.encoding:
            logo = """
╔╦╗╦═╗╦╔═╗ ╔═╗  ┬ ┬┌─┐
 ║ ╠╦╝║║═╬╗╚═╗  ├─┤├┤ 
 ╩ ╩╚═╩╚═╝╚╚═╝  ┴ ┴└  
TRIQS: Hartree-Fock solver
"""
        else:
            logo = """
 _____ ___ ___ ___  ___   _  _ ___
|_   _| _ \_ _/ _ \/ __| | || | __|
  | | |   /| | (_) \__ \ | __ | _|
  |_| |_|_\___\__\_\___/ |_||_|_|

    TRIQS: Hartree-Fock solver
"""
        return logo
    
# ----------------------------------------------------------------------
class HartreeSolver(HartreeFockSolver):

    """ Hartree solver for local interactions.

    Parameters
    ----------

    e_k : TRIQS Greens function on a Brillouin zone mesh
        Single-particle dispersion.

    beta : float
        inverse temperature

    H_int : TRIQS Operator instance
        Local interaction Hamiltonian

    gf_struct : list of lists of block index and list of internal orbital indices
        gf_struct fixing orbital order between e_k and H_int

    mu0 : float
        chemical potential

    mu_min, mu_max : float, optional
        range for chemical potential search.

    """

    # ------------------------------------------------------------------
    def mat2vec(self, mat):
        assert( len(mat.shape) == 2 )
        return np.diag(mat).real

    # ------------------------------------------------------------------
    def vec2mat(self, vec):
        assert( len(vec.shape) == 1 )
        return np.diag(vec.real)

    # ------------------------------------------------------------------
    def logo(self):
        if 'UTF' in sys.stdout.encoding:
            logo = """
╔╦╗╦═╗╦╔═╗ ╔═╗  ┬ ┬
 ║ ╠╦╝║║═╬╗╚═╗  ├─┤
 ╩ ╩╚═╩╚═╝╚╚═╝  ┴ ┴
TRIQS: Hartree solver
"""
        else:
            logo = """
 _____ ___ ___ ___  ___   _  _
|_   _| _ \_ _/ _ \/ __| | || |
  | | |   /| | (_) \__ \ | __ |
  |_| |_|_\___\__\_\___/ |_||_|

    TRIQS: Hartree solver
"""
        return logo

# ----------------------------------------------------------------------
