
import itertools
import numpy as np

# ----------------------------------------------------------------------    

from triqs.operators import Operator
from triqs.operators.util.U_matrix import U_matrix_kanamori
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs.operators.util.op_struct import set_operator_structure

# ----------------------------------------------------------------------    

from triqs_tprf.OperatorUtils import quartic_permutation_symmetrize
from triqs_tprf.OperatorUtils import fundamental_operators_from_gf_struct
from triqs_tprf.OperatorUtils import get_quadratic_operator,  quadratic_matrix_from_operator
from triqs_tprf.OperatorUtils import quartic_tensor_from_operator, operator_from_quartic_tensor
from triqs_tprf.OperatorUtils import operator_single_particle_transform, relabel_operators

# ----------------------------------------------------------------------    

from triqs_tprf.rpa_tensor import get_rpa_tensor
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors
from triqs_tprf.rpa_tensor import split_quartic_tensor_in_charge_and_spin
from triqs_tprf.rpa_tensor import quartic_tensor_from_charge_and_spin

# ----------------------------------------------------------------------    
def print_tensors(T1, T2):

    assert( T1.shape == T2.shape )
    for idxs in itertools.product(*[ list(range(x)) for x in T1.shape ]):
        print(idxs, T1[idxs], T2[idxs])
        
# ----------------------------------------------------------------------    
if __name__ == '__main__':

    U = 1.0
    J = 0.3
    
    orb_names = ['xy', 'xz', 'yz']
    spin_names = ['up', 'do']
    norb = len(orb_names)

    gf_struct = set_operator_structure(spin_names, orb_names, True) # orbital off diag
    fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)
    
    U_ab, UPrime_ab = U_matrix_kanamori(
        n_orb=len(orb_names), U_int=U, J_hund=J)
    
    H_int = h_int_kanamori(
        spin_names, orb_names, U_ab, UPrime_ab, J_hund=J,
        off_diag=True, map_operator_structure=None, H_dump=None) # orgital offdiag

    U_abcd = get_rpa_tensor(H_int, fundamental_operators)
    
    U_c, U_s = split_quartic_tensor_in_charge_and_spin(U_abcd)

    U_abcd_ref = quartic_tensor_from_charge_and_spin(U_c, U_s)

    U_c_ref, U_s_ref = kanamori_charge_and_spin_quartic_interaction_tensors(
        norb, U, U - 2*J, J, J)

    np.testing.assert_array_almost_equal(U_abcd, U_abcd_ref)
    np.testing.assert_array_almost_equal(U_c, U_c_ref)
    np.testing.assert_array_almost_equal(U_s, U_s_ref)
    
