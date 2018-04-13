# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import Gf
from pytriqs.gf import Idx
from pytriqs.gf import MeshImFreq
from pytriqs.gf import MeshProduct

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

# ----------------------------------------------------------------------
def chi_w0r_from_chi_tr_np_trapz(chi_tr):
    
    tmesh = chi_tr.mesh.components[0]
    rmesh = chi_tr.mesh.components[1]

    beta = tmesh.beta
    tau = np.array([float(t) for t in tmesh])
    
    wmesh = MeshImFreq(beta=beta, S='Boson', n_max=1)

    #print 'tau =', tau
    chi00_wr = Gf(mesh=MeshProduct(wmesh, rmesh), target_shape=chi_tr.target_shape)
    chi00_wr[Idx(0), :].data[:] = np.trapz(chi_tr.data, x=tau, axis=0)

    return chi00_wr

# ----------------------------------------------------------------------
def cluster_mesh_fourier_interpolation(k, chiwr):

    assert( len(k.shape) == 2 )
    assert( k.shape[1] == 3 )
    
    rmesh = chiwr.mesh.components[1]
    r = np.array([r.value for r in rmesh])
    
    exp_fact = np.exp(np.sum(-1.j * k[:, None, :]  * r[None, :, :], axis=-1))

    chi00wk_data = np.einsum('wrabcd,kr->wkabcd', chiwr.data, exp_fact)

    return chi00wk_data
