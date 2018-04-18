# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import NearestNDInterpolator        

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
    k_vec_rel = get_relative_k_from_absolute(k_vec, bz.units())
    
    ek.data[:] = tb_lattice.hopping(k_vec_rel.T).transpose(2, 0, 1)

    return ek

# ----------------------------------------------------------------------
def get_relative_k_from_absolute(k_vec, units):
    k_vec_rel = np.dot(np.linalg.inv(units).T, k_vec.T).T
    return k_vec_rel

# ----------------------------------------------------------------------
def get_kidx_from_k_vec_relative(k_vec_rel, periodization_matrix):
    nk = np.diag(periodization_matrix)
    kidx = np.array(np.round(k_vec_rel * nk[None, :]), dtype=np.int)
    return kidx

# ----------------------------------------------------------------------
def get_k_components_from_k_vec(k_vec, periodization_matrix):

    dim = 3
    nk = np.diag(periodization_matrix)
    #print 'nk =', nk

    shape_4 = [dim] + list(nk)
    shape_2 = [dim] + [np.prod(nk)]

    k_vec = k_vec.swapaxes(0, -1)    
    k_vec = k_vec.reshape(shape_4)

    # -- cut out values for each axis
    
    k_out = []
    for axis in xrange(dim):
        cut = [0]*dim
        cut[axis] = slice(None)
        cut = [axis] + cut
        k_out.append(k_vec[tuple(cut)])
        
    return tuple(k_out)

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

# ----------------------------------------------------------------------
def get_abs_k_chi_interpolator(chi, w, bz, extend_bz=[0]):

    k_mat = bz.units()
    bzmesh = chi.mesh.components[1]

    k_vec = np.array([k.value for k in bzmesh])
    values = np.squeeze(chi[w, :].data.real)
    
    # -- Extend with points beyond the first bz

    k_vec_ext = []
    values_ext = []
    for k_shift in itertools.product(extend_bz, repeat=3):
        k_shift = np.dot(k_mat.T, k_shift)        
        k_vec_ext.append( k_vec + k_shift[None, :] )
        values_ext.append(values)

    k_vec = np.vstack(k_vec_ext)
    values = np.hstack(values_ext)

    interp = LinearNDInterpolator(k_vec, values, fill_value=float('nan'))
    
    return interp
    
# ----------------------------------------------------------------------
def get_rel_k_chi_interpolator(chi, w, bz, periodization_matrix,
                               extend_boundary=True, interpolator='regular'):

    k_mat = bz.units()
    bzmesh = chi.mesh.components[1]
    nk = np.diag(periodization_matrix)

    k_vec = np.array([k.value for k in bzmesh])
    values = np.squeeze(chi[w, :].data.real)

    k_vec_rel = get_relative_k_from_absolute(k_vec, bz.units())
    k_idx = get_kidx_from_k_vec_relative(k_vec_rel, periodization_matrix)

    kx, ky, kz = get_k_components_from_k_vec(k_vec_rel, periodization_matrix)

    if extend_boundary:
        values, k_vec_rel, (kx, ky, kz) = \
            extend_data_on_boundary(values, periodization_matrix)
    else:
        values = values.reshape(nk)

    # -- select interpolator type
        
    if interpolator is 'regular':
        interp = RegularGridInterpolator(
            (kx, ky, kz), values, fill_value=float('nan'), bounds_error=False)
    elif interpolator is 'nearest':
        interp = NearestNDInterpolator(k_vec_rel, values.flatten())
    elif interpolator is 'linear':
        interp = LinearNDInterpolator(k_vec_rel, values.flatten(), fill_value=float('nan'))
    else:
        raise NotImplementedError
        
    return interp

# ----------------------------------------------------------------------
def extend_data_on_boundary(values, periodization_matrix):

    # -- extended cube

    nk = np.diag(periodization_matrix)

    nk_ext = nk + 1

    coords = [ np.arange(0, n) for n in nk_ext ]
    Coords = np.meshgrid(*coords, indexing='ij')
    Coords_mod = [ np.mod(x, n) for x, n in zip(Coords, nk) ]

    values_ext = values.reshape(nk)[tuple([ X.flatten() for X in Coords_mod])]
    values_ext = values_ext.reshape(nk_ext)

    # -- compute kidx_ext

    k_idx_ext = np.array([ X.flatten() for X in Coords ]).T
    k_vec_rel_ext = np.array(k_idx_ext, dtype=np.float) / nk[None, :]
    kxe, kye, kze = get_k_components_from_k_vec(k_vec_rel_ext, np.diag(nk_ext))

    return values_ext, k_vec_rel_ext, (kxe, kye, kze)
