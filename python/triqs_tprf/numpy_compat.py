
import numpy as np
from packaging.version import parse

def is_numpy_newer_than(version):
    return parse(np.__version__) > parse(version)

def np_linalg_func(A, func, version='1.8.0'):
    if is_numpy_newer_than(version): Ainv = func(A)
    else: Ainv = np.array(map(func, A))
    return Ainv

def np_inv(A):
    return np_linalg_func(A, np.linalg.inv)

def np_eigh(A):
    return np_linalg_func(A, np.linalg.eigh)

def np_eigvalsh(A):
    return np_linalg_func(A, np.linalg.eigvalsh)
