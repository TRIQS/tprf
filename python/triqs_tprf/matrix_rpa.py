# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019, S. Käser
# Author: S. Käser
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

# ----------------------------------------------------------------------

""" This set of functions implements the matrix RPA as preseneted in mutliple papers.
(Note that there is still the contradiction, that we use chi*U indest of U*chi.)
Four rank tensors in c^+ccc^+ (susceptibilites) or cc^+c^+c (vertices) order
are reshaped to matrices, and the RPA equations are simple matrix multiplications.
This is used to test if TPRF and matrix RPA yields the same results.
"""

# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------
def tprf_order_to_matrix_rpa_order(tensor):
    """Bring the order used in tprf, i.e. c^+ccc^+ for susceptibilites
    and cc^+c^+c for interactions into the order of matrix RPA, i.e.
    c^+ccc^+ for susceptibilites and cc^+c^+c for interacions.
    """

    return tensor.swapaxes(-1,-2)

# ----------------------------------------------------------------------
def lose_spin_degree_of_freedom(gf, rank=2, spin_fast=True):
    """Only keep the up/down spin elements of a tensor

    Parameters:

    gf: np.ndarray, the last rank dimenions must the orbitals
    rank: int, rank of the tensor
    spin_fast: bool, True if spin is the fast index, e.g.
                    xz up, xz down, xy up, xy down, yz up, yz down,
                or False if spin is the slow index, e.g.
                    xz up, xy up, yz up, xz down, xy down, yz down.
    """

    norb = gf.shape[-1] // 2

    idx = (len(gf.shape)-rank) * [slice(None)]
    
    if spin_fast:
        idx.extend(rank*[slice(None, None, 2)])

    else:
        idx.extend(rank*[slice(norb, None)])
    
    spin_independent_gf = gf[tuple(idx)]

    return spin_independent_gf

# ----------------------------------------------------------------------
def get_rpa_us_tensor(norb, U, Up, J, Jp):

    """Returns interaction tensor for spin susceptibility 
    cc^+c^+c order

    Parameters:

    norb : int, number of orbitals of the model
    U : float, intra orbital interaction
    Up : float, inter orbital interaction
    J : float, Hund's coupling
    Jp : float, pair hopping and spin flip
    """

    U_abcd = np.zeros(4*[norb], dtype=complex)

    for a,b,c,d in itertools.product(list(range(norb)), repeat=4):

        if a == b == c == d:
            U_abcd[a,b,c,d] = U

        elif a == c != b == d:
            U_abcd[a,b,c,d] = Up

        elif a == b != c == d:
            U_abcd[a,b,c,d] = J

        elif a == d != b == c:
            U_abcd[a,b,c,d] = Jp

    return U_abcd
    
# ----------------------------------------------------------------------
def get_rpa_uc_tensor(norb, U, Up, J, Jp):

    """Returns interaction tensor for charge susceptibility 
    cc^+c^+c order

    Parameters:

    norb : int, number of orbitals of the model
    U : float, intra orbital interaction
    Up : float, inter orbital interaction
    J : float, Hund's coupling
    Jp : float, pair hopping and spin flip
    """

    U_abcd = np.zeros(4*[norb], dtype=complex)

    for a,b,c,d in itertools.product(list(range(norb)), repeat=4):

        if a == b == c == d:
            U_abcd[a,b,c,d] = U

        elif a == c != b == d:
            U_abcd[a,b,c,d] = -Up + 2*J

        elif a == b != c == d:
            U_abcd[a,b,c,d] = 2*Up - J

        elif a == d != b == c:
            U_abcd[a,b,c,d] = Jp

    return U_abcd

# ----------------------------------------------------------------------
def tensor_to_matrix(tensor):

    """Tranforms a 4-rank tensor into a matrix using the matrix-RPA notation 

    Parameters:

    tensor: np.ndarray, last 4 dimensions must be the orbitals
    """

    norb = tensor.shape[-1]
    non_orb_idx = tensor.shape[:-4]
    
    idx = list(range(norb**2))
    diag = idx[::norb+1]
    idx = diag + list(set(idx) - set(diag))
    
    left = [slice(None)]*len(non_orb_idx) + [np.array(idx)]
    right = [slice(None)]*(len(non_orb_idx)+1) + [np.array(idx)]
    
    matrix = tensor.reshape(non_orb_idx+(norb**2,)*2)[tuple(left)][tuple(right)]

    return matrix

# ----------------------------------------------------------------------
def matrix_to_tensor(matrix):

    """Tranforms a matrix-RPA notation matrix into a 4-rank tensor

    Parameters:

    matrix: np.ndarray, last 2 dimensions must be the orbitals
    """

    norb = int(np.sqrt(matrix.shape[-1]))
    non_orb_idx = matrix.shape[:-2]

    idx = list(range(norb**2))
    diag = idx[::norb+1]
    idx = diag + list(set(idx) - set(diag))
    idx = [idx.index(l) for l in range(len(idx))]
    
    left = [slice(None)]*len(non_orb_idx) + [np.array(idx)]
    right = [slice(None)]*(len(non_orb_idx)+1) + [np.array(idx)]
    
    tensor = matrix[tuple(left)][tuple(right)].reshape(non_orb_idx + (norb,)*4)

    return tensor 

# ----------------------------------------------------------------------
def matrix_rpa(chi0_matrix, u_matrix):

    """Calculates the matrix RPA equation

    Parameters:

    chi0: np.ndarray, matrix RPA notation with the last two dimensions as orbitals
    u: np.ndarray, matrix RPA notation with the last two dimensions as orbitals
    """

    norb = int(np.sqrt(chi0_matrix.shape[-1]))
    
    u_chi0 = np.einsum('...ij,...jk->...ik', chi0_matrix, u_matrix)

    inverse_matrix = np.linalg.inv(np.eye(norb**2) - u_chi0)
    
    chi_rpa_matrix = np.einsum('...ij,...jk->...ik', inverse_matrix, chi0_matrix)

    return chi_rpa_matrix
    
# ----------------------------------------------------------------------
def chi_rpa_spin(chi0, us):

    """Calculates the spin susceptibility tensor

    Parameters:

    chi0: np.ndarray, 4-rank tensor with the last 4 dimensions as orbitals,
            ordering c^+ccc^+
    us: np.ndarray, 4-rank tensor with the last 4 dimensions as orbitals
            ordering cc^+c^+c
    """
    
    chi0_matrix = tensor_to_matrix(chi0)
    us_matrix = tensor_to_matrix(us)

    chi_rpa_spin_matrix = matrix_rpa(chi0_matrix, us_matrix)

    chi_rpa_spin_tensor = matrix_to_tensor(chi_rpa_spin_matrix)

    return chi_rpa_spin_tensor

# ----------------------------------------------------------------------
def chi_rpa_charge(chi0, uc):

    """Calculates the spin susceptibility tensor

    Parameters:

    chi0: np.ndarray, 4-rank tensor with the last 4 dimensions as orbitals,
            ordering c^+ccc^+
    uc: np.ndarray, 4-rank tensor with the last 4 dimensions as orbitals
            ordering cc^+c^+c
    """
    
    chi0_matrix = tensor_to_matrix(chi0)
    uc_matrix = tensor_to_matrix(uc)
    
    # -- Minus infront of uc to get the right equation for the charge susceptibility
    chi_rpa_charge_matrix = matrix_rpa(chi0_matrix, -uc_matrix)

    chi_rpa_charge_tensor = matrix_to_tensor(chi_rpa_charge_matrix)

    return chi_rpa_charge_tensor
