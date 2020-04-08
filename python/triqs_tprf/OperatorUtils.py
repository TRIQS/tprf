
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

""" Utilities for working with Triqs second quantized operator
expressions.

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import itertools
import numpy as np

from pytriqs.operators import c, c_dag, Operator, dagger

# ----------------------------------------------------------------------    
def fundamental_operators_from_gf_struct(gf_struct):

    """ Return a list of annihilation operators with the quantum numbers 
    defined in the gf_struct """
    
    fundamental_operators = []
    for block_name, indices in gf_struct:
        for index in indices:
            fundamental_operators.append(
                c(block_name, index) )

    return fundamental_operators

# ----------------------------------------------------------------------
def get_quadratic_operator(h, fundamental_operators):

    # -- Check Hermicity
    np.testing.assert_array_almost_equal(h, h.T.conj())
    
    H = Operator(0.)
    for idx1, o1 in enumerate(fundamental_operators):
        o1 = dagger(o1)
        for idx2, o2 in enumerate(fundamental_operators):
            H += h[idx1, idx2] * o1 * o2

    return H
            
# ----------------------------------------------------------------------
def op_is_fundamental(op):
    
    """ Check for single "fundamental" operator in operator expression. """

    monomial_sum = list(op)    
    if len(monomial_sum) != 1:
        return False

    monomial = monomial_sum[0]
    operator_list, prefactor = monomial

    if len(operator_list) != 1:
        return False

    return True

# ----------------------------------------------------------------------
def op_serialize_fundamental(op):
    
    """ Return a tuple specifying dagger and operator indices of a 
    fundamental (single) creation/annihilation operator. """

    assert( op_is_fundamental(op) )
    dag, idxs = list(op)[0][0][0]
    return dag, tuple(idxs)
    
# ----------------------------------------------------------------------
def get_operator_index_map(fundamental_operators, include_dag=False):

    """ Construct list of tuples of orbital index pairs from list of
    fundamental operators. """

    op_idx_map = []
    for op in fundamental_operators:
        dag, idxs = op_serialize_fundamental(op)
        if include_dag:
            op_idx_map.append((dag, tuple(idxs)))
        else:
            assert( dag is False ) # Only accept annihilation ops
            op_idx_map.append(tuple(idxs))

    return op_idx_map

# ----------------------------------------------------------------------
def quadratic_matrix_from_operator(op, fundamental_operators):

    """ Build matrix of quadratic operator terms spanned by the list
    of fundamental_operators in given operator expression op. """

    op_idx_map = get_operator_index_map(fundamental_operators)
    op_idx_set = set(op_idx_map)

    nop = len(fundamental_operators)
    h_quad = np.zeros((nop, nop), dtype=np.complex)
    
    for term in op:
        op_list, prefactor = term
        if len(op_list) == 2:
            (d1, t1), (d2, t2) = list(op_list)
            t1, t2 = tuple(t1), tuple(t2)
            assert( d1 is True and d2 is False ) # no anomalous terms, yet...
            if t1 in op_idx_set and t2 in op_idx_set:
                i, j = op_idx_map.index(t1), op_idx_map.index(t2)
                h_quad[i, j] = prefactor

    return h_quad

# ----------------------------------------------------------------------
def quartic_permutation_symmetrize(U):

    U = 0.5 * ( U - np.swapaxes(U, 0, 1) )
    U = 0.5 * ( U - np.swapaxes(U, 2, 3) )

    return U

# ----------------------------------------------------------------------
def quartic_conjugation_symmetrize(U):

    U = 0.5 * ( U + np.swapaxes(np.swapaxes(U, 1, 2), 0, 3).conj() )

    return U

# ----------------------------------------------------------------------
def quartic_pauli_symmetrize(U):

    N = U.shape[0]
    assert( U.shape == tuple([N]*4) )

    for n in range(N):
        U[n, n, :, :] = 0
        U[:, :, n, n] = 0

    return U

# ----------------------------------------------------------------------
def symmetrize_quartic_tensor(U, conjugation=False):

    N = U.shape[0]
    assert( U.shape == tuple([N]*4) )
    
    # -- Permutation symmetries

    #U = 0.5 * ( U - np.swapaxes(U, 0, 1) )
    #U = 0.5 * ( U - np.swapaxes(U, 2, 3) )

    U = quartic_permutation_symmetrize(U)

    # -- Conjugation symmetry

    if conjugation:
        #U = 0.5 * ( U + np.swapaxes(np.swapaxes(U, 1, 2), 0, 3).conj() )
        U = quartic_conjugation_symmetrize(U)

    # -- Pauli principle
    
    #for n in xrange(N):
    #    U[n, n, :, :] = 0
    #    U[:, :, n, n] = 0

    U = quartic_pauli_symmetrize(U)
    
    return U
    
# ----------------------------------------------------------------------
def quartic_tensor_from_operator(op, fundamental_operators,
                                 perm_sym=False):

    # -- Convert fundamental operators back and forth from index

    op_idx_map = get_operator_index_map(fundamental_operators)
    op_idx_set = set(op_idx_map)

    nop = len(fundamental_operators)
    h_quart = np.zeros((nop, nop, nop, nop), dtype=np.complex)
    
    for term in op:
        op_list, prefactor = term
        if len(op_list) == 4:

            d, t = list(zip(*op_list)) # split in two lists with daggers and tuples resp
            t = [tuple(x) for x in t]

            # check creation/annihilation order
            assert( d == (True, True, False, False) ) 
            
            if all([ x in op_idx_set for x in t ]):
                i, j, k, l = [ op_idx_map.index(x) for x in t ]

                if perm_sym:
                    # -- all pair wise permutations:
                    # -- i <-> j and k <-> l
                    # -- with fermionic permutation signs

                    h_quart[i, j, k, l] = +0.25 * prefactor
                    h_quart[j, i, k, l] = -0.25 * prefactor
                    h_quart[j, i, l, k] = +0.25 * prefactor
                    h_quart[i, j, l, k] = -0.25 * prefactor
 
                else:

                    h_quart[i, j, k, l] = prefactor

    return h_quart

# ----------------------------------------------------------------------
def operator_from_quartic_tensor(h_quart, fundamental_operators):

    # -- Convert fundamental operators back and forth from index

    op_idx_map = get_operator_index_map(fundamental_operators)
    op_idx_set = set(op_idx_map)

    nop = len(fundamental_operators)

    H = Operator(0.)

    for t in itertools.product(enumerate(fundamental_operators), repeat=4):
        idx, ops = list(zip(*t))
        o1, o2, o3, o4 = ops
        o1, o2 = dagger(o1), dagger(o2)

        H += h_quart[idx] * o1 * o2 * o3 * o4

    return H    

# ----------------------------------------------------------------------
def operator_single_particle_transform(op, U, fundamental_operators):
    
    """ Transform the operator op according to the single particle 
    transform matrix U defined in the subspace of operators listed in 
    fundamental_operators. """

    # -- Convert fundamental operators back and forth from index

    op_idx_map = get_operator_index_map(fundamental_operators)
    op_idx_set = set(op_idx_map)

    # -- Transformed creation operator

    def c_transf(s, i):
        if (s, i) not in op_idx_set:
            return c(s, i)

        k = op_idx_map.index((s, i))
        
        ret = Operator()
        for l in range(U.shape[0]):
            op_idx = op_idx_map[l]
            ret += U[k, l] * c(*op_idx)

        return ret

    # -- Precompute transformed operators
    
    op_trans_dict = {}
    for fop in fundamental_operators:
        dag, idxs = op_serialize_fundamental(fop)
        op_trans_dict[(dag, idxs)] = c_transf(*idxs)
        op_trans_dict[(not dag, idxs)] = dagger(c_transf(*idxs))
            
    # -- Loop over given operator and substitute operators
    # -- fundamental_operators with the transformed operators
    
    op_trans = Operator()
    
    for term in op:
        op_factor = Operator(1.)
        for factor in term:
            if type(factor) is list:
                for dag, idxs in factor:
                    tup = (dag, tuple(idxs))
                    if tup in list(op_trans_dict.keys()):
                        op_factor *= op_trans_dict[tup]
                    else:
                        op_factor *= {False:c, True:c_dag}[dag](*idxs)
                    
            else: # constant prefactor
                op_factor *= factor

        op_trans += op_factor

    return op_trans

# ----------------------------------------------------------------------
def relabel_operators(op, from_list, to_list):
    
    """ Transform the operator op according to the single particle 
    transform matrix U defined in the subspace of operators listed in 
    fundamental_operators. """

    op_idx_map = get_operator_index_map(from_list, include_dag=True)
    op_idx_set = set(op_idx_map)

    from_list, to_list = list(from_list), list(to_list)

    # -- Loop over given operator and substitute operators
    # -- in from_list to to_list
    
    op_trans = Operator()
    
    for term in op:
        op_factor = Operator(1.)
        for factor in term:
            if type(factor) is list:
                for dag, idxs in factor:
                    tup = (dag, tuple(idxs))
                    if tup in op_idx_set:
                        op_factor *= to_list[op_idx_map.index(tup)]
                    else:
                        op_factor *= {False:c, True:c_dag}[dag](*idxs)
                    
            else: # constant prefactor
                op_factor *= factor

        op_trans += op_factor

    return op_trans

# ----------------------------------------------------------------------
def is_operator_composed_of_only_fundamental_operators(
        op, fundamental_operators):

    """ Return `True` if the operator `op` only contains operators belonging
    to the fundamental operator set `fundamental_operators`, else return
    `False`."""
    
    is_fundamental = True

    for term in op:
        op_list, prefactor = term
        d, t = list(zip(*op_list)) # split in two lists with daggers and tuples resp
        t = [tuple(x) for x in t]
        for bidx, idx in t:
            if c(bidx, idx) not in fundamental_operators:
                is_fundamental = False

    return is_fundamental
