# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import Gf

# ----------------------------------------------------------------------
def ek_tb_dispersion_on_bzmesh(tb_lattice, bzmesh, bz):

    """ Evaluate dispersion on bzmesh from tight binding model. """

    n_orb = tb_lattice.NOrbitalsInUnitCell
    ek = Gf(mesh=bzmesh, target_shape=[n_orb]*2)

    k_vec = np.array([k.value for k in bzmesh])

    k_mat = bz.units()
    k_vec = np.dot(np.linalg.inv(k_mat).T, k_vec.T).T

    ek.data[:] = tb_lattice.hopping(k_vec.T).transpose(2, 0, 1)

    return ek
