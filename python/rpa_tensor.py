
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from pyed.OperatorUtils import quartic_tensor_from_operator
from pyed.OperatorUtils import operator_from_quartic_tensor

from pyed.OperatorUtils import symmetrize_quartic_tensor

from pyed.OperatorUtils import quartic_permutation_symmetrize
from pyed.OperatorUtils import quartic_conjugation_symmetrize
from pyed.OperatorUtils import quartic_pauli_symmetrize

# ----------------------------------------------------------------------
def fundamental_operators_from_gf_struct(gf_struct):

    fundamental_operators = []
    for block, idxs in gf_struct:
        for idx in idxs:
            fundamental_operators.append( c(block, idx) )

    return fundamental_operators

# ----------------------------------------------------------------------
def get_rpa_tensor(H_int, fundamental_operators):

    U_int = quartic_tensor_from_operator(H_int, fundamental_operators)
    U_int = 4 * quartic_permutation_symmetrize(U_int)

    # -- Group in Gamma order cc^+cc^+ ( from c^+c^+cc )
    
    U_abcd = np.zeros_like(U_int)
    for a, b, c, d in itertools.product(range(U_int.shape[0]), repeat=4):
        U_abcd[a, b, c, d] = U_int[b, d, a, c]

    return U_abcd
