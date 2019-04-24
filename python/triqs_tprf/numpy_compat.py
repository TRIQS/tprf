
""" Helper functions that enable the numpy >= 1.8.0
behavior of np.linalg functions.

If a three dimensional array are passed to, inv, eigh, eigvalsh, etc.
the linear algebra operation is applied to the last two indices. 

These routines use this behavior if the numpy version is high enough
and uses a fallback if that is not the case.

Author: H. U.R. Strand (2019)
"""

import numpy as np
#from packaging.version import parse
from distutils.version import StrictVersion as parse

def is_numpy_newer_than(version):
    return parse(np.__version__) > parse(version)

def np_linalg_func(arr, func, version='1.8.0'):
    if is_numpy_newer_than(version): arr_out = func(arr)
    else: arr_out = np.array(map(func, arr))
    return arr_out

def np_inv(A):
    return np_linalg_func(A, np.linalg.inv)

def np_eigh(A):
    if is_numpy_newer_than(version):
        E, U = np.linalg.eigh(A)
    else:
        N, M, K = A.shape
        assert(M == K)
        E = np.empty((N, M))
        U = np.empty((N, M, M))        
        for i in xrange(N):
            E[i], U[i] = np.linalg.eigh(A[i])
            
    return E, U

def np_eigvalsh(A):
    return np_linalg_func(A, np.linalg.eigvalsh)
