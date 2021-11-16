# ----------------------------------------------------------------------

""" Symmetrize a randomly filled Green's function in frequency, momentum,
    and orbital and test if it was done proberly.

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import Gf, MeshImFreq, MeshBrZone, MeshProduct
from triqs.lattice import BrillouinZone, BravaisLattice
from triqs_tprf.ParameterCollection import *

# ----------------------------------------------------------------------

from triqs_tprf.symmetries import enforce_symmetry, check_symmetry, _invert_momentum

# ----------------------------------------------------------------------

p = ParameterCollection(beta = 10,
                        nw = 10,
                        nk = 4,
                        norb = 2,)

wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nw)

cell = np.eye(3)
bl = BravaisLattice(cell)
bz = BrillouinZone(bl)
kmesh = MeshBrZone(bz, p.nk * np.eye(3, dtype=int))

gf = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=2*(p.norb,))
gf.data[:] = np.random.rand(*gf.data.shape)

# -- Eexception handling
try:
    enforce_symmetry(gf, "something", "odd")
except ValueError as error:
    if not str(error) == "No symmetrize function for this variable exists.":
        raise Exception("Wrong exception was raised: \n %s"%error)
else:
    raise Exception("Function call should have failed.")

try:
    enforce_symmetry(gf, "frequency", "weird")
except ValueError as error:
    if not str(error) == "Symmetry can only be 'even' or 'odd'.":
        raise Exception("Wrong exception was raised: \n %s"%error)
else:
    raise Exception("Function call should have failed.")

# -- Frequency
even_freq_gf = enforce_symmetry(gf, "frequency", "even")
np.testing.assert_equal(+1, check_symmetry(even_freq_gf)['frequency'])

odd_freq_gf = enforce_symmetry(gf, "frequency", "odd")
np.testing.assert_equal(-1, check_symmetry(odd_freq_gf)['frequency'])

# -- Momentum 
momentum_mesh = [16, 32, 5]
momenta = [ [0, 0, 0],
            [1, 4, 0],
            [8, 16, 3],
          ]
inv_momenta = [ [0, 0, 0],
                [15, 28, 0],
                [8, 16, 2],
              ]
for momentum, inv_momentum in zip(momenta, inv_momenta):
    inv_momentum_test = _invert_momentum(momentum, momentum_mesh)
    if not np.equal(inv_momentum, inv_momentum_test).all:
        raise Exception("The function '_invert_momentum' does not behave as"
                        " expected.")

even_momentum_gf = enforce_symmetry(gf, "momentum", "even")
np.testing.assert_equal(+1, check_symmetry(even_momentum_gf)["momentum"])

odd_momentum_gf = enforce_symmetry(gf, "momentum", "odd")
np.testing.assert_equal(-1, check_symmetry(odd_momentum_gf)["momentum"])

# -- Orbital

even_orbital_gf = enforce_symmetry(gf, "orbital", "even")
np.testing.assert_equal(+1, check_symmetry(even_orbital_gf)["orbital"])

odd_orbital_gf = enforce_symmetry(gf, "orbital", "odd")
np.testing.assert_equal(-1, check_symmetry(odd_orbital_gf)["orbital"])

# -- Combination
variables = ["frequency", "momentum", "orbital"]
avail_symmetries = ["even", "odd", None]
for symmetries in itertools.product(avail_symmetries, repeat=len(variables)):

    # Remove the corresponding variable to None that it does not get symmetrized
    variables_copy = list(variables)
    symmetries_copy = list(symmetries)
    while None in symmetries_copy:
        idx = symmetries_copy.index(None)
        del symmetries_copy[idx]
        del variables_copy[idx]

    symmetrized_gf = enforce_symmetry(gf, variables_copy, symmetries_copy)
    
    translate_symmetries = {"even" : +1, "odd" : -1, None : None}
    expected_symmetries = {variable : translate_symmetries[symmetry] \
                    for (variable, symmetry) in zip(variables, symmetries)}
        
    produced_symmetries = check_symmetry(symmetrized_gf)

    if not expected_symmetries == produced_symmetries:
        raise AssertionError("Incorrect symmetries were produced")

print("It's all good honey")
