
""" Helper functions that enable the numpy >= 1.8.0
behavior of np.linalg functions.

If a three dimensional array are passed to, inv, eigh, eigvalsh, etc.
the linear algebra operation is applied to the last two indices. 

These routines use this behavior if the numpy version is high enough
and uses a fallback if that is not the case.

Author: H. U.R. Strand (2019)
"""

import itertools
import numpy as np
from packaging.version import parse

def is_numpy_newer_than(version):
    return parse(np.__version__) > parse(version)

def np_linalg_func(arr, func, version='1.8.0'):
    #if False:
    if is_numpy_newer_than(version):
        arr_out = func(arr)
    else:
        arr_out = np.empty_like(arr)
        ranges = [range(N) for N in arr.shape[:-2]]
        for idx in itertools.product(*ranges):
            arr_out[idx] = func(arr[idx])

    #arr_out_ref = func(arr)
    #np.testing.assert_array_almost_equal(arr_out, arr_out_ref)
    return arr_out

def np_inv(A):
    return np_linalg_func(A, np.linalg.inv)

def np_eigh(A):
    #if False:
    if is_numpy_newer_than('1.8.0'):
        E, U = np.linalg.eigh(A)
    else:
        N, M, K = A.shape
        assert(M == K)
        E = np.empty((N, M))
        U = np.empty((N, M, M), dtype=complex) 
        for i in range(N):
            E[i], U[i] = np.linalg.eigh(A[i])
            
    #E_ref, U_ref = np.linalg.eigh(A)
    #np.testing.assert_array_almost_equal(E, E_ref)
    #np.testing.assert_array_almost_equal(U, U_ref)
    return E, U

def np_eigvalsh(arr):
    #if False:
    if is_numpy_newer_than('1.8.0'):
        arr_out = np.linalg.eigvalsh(arr)
    else:
        arr_out = np.empty_like(np.squeeze(arr[...,0]))
        ranges = [range(N) for N in arr.shape[:-2]]
        for idx in itertools.product(*ranges):
            arr_out[idx] = np.linalg.eigvalsh(arr[idx])

    #arr_out_ref = np.linalg.eigvalsh(arr)
    #np.testing.assert_array_almost_equal(arr_out, arr_out_ref)
    return arr_out
