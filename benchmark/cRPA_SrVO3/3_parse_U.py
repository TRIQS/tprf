"""
Read the local interaction tensor from UIJKL
and compare to the momentum dependent interaction UIJKL_Q after summing over Q.

Author: H. U.R. Strand
"""

import itertools
import numpy as np

from eval_U import read_uijkl, red_to_2ind
from vasp_crpa_parsers import read_vasp_crpa_momentum_space_interaction_to_ndarray

print '='*72
print '--> Reading UIJKL file'
print '='*72

U_ijkl = read_uijkl('./crpa/UIJKL')
U1, U2, U3, U4 = red_to_2ind(U_ijkl, verbose=True)

# -- Read q-dependent interaction

print '='*72
print '--> Reading UIJKL_Q and VIJKL_Q files'
print '='*72

UR_Q, Q = read_vasp_crpa_momentum_space_interaction_to_ndarray('./crpa', 'UIJKL_Q_full.q*')
VR_Q, Q = read_vasp_crpa_momentum_space_interaction_to_ndarray('./crpa', 'VIJKL_Q_full.q*')
VRR_Q, Q = read_vasp_crpa_momentum_space_interaction_to_ndarray('./crpa', 'VIJKL_Q_redu.q*')
U_Q = UR_Q + ( VR_Q - VRR_Q )

U_ijkl_ref = np.sum(U_Q, axis=0) / Q.shape[0]
U1_ref, U2_ref, U3_ref, U4_ref = red_to_2ind(U_ijkl_ref, verbose=True)

print '='*72
print '--> Full tensor comparison'
print '='*72

for i, j, k, l in itertools.product(range(3), repeat=4):
    a, b = U_ijkl[i, j, k, l].real, U_ijkl_ref[i, j, k, l].real
    c, d = U_ijkl[i, j, k, l].imag, U_ijkl_ref[i, j, k, l].imag
    print i, j, k, l, 'RE %+6.6f %+6.6f (diff %+6.6f) -- IM  %+6.6f %+6.6f (diff %+6.6f)' % (a, b, a - b, c, d, c - d)

diff = np.max(np.abs( U_ijkl - U_ijkl_ref ))
print 'max(abs(diff)) =', diff
