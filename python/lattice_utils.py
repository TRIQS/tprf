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
from pytriqs.gf import MeshBrillouinZone

from pytriqs.lattice import BrillouinZone

# ----------------------------------------------------------------------

from triqs_tprf.logo import tprf_banner

from triqs_tprf.lattice import chi00_wk_from_ek

from triqs_tprf.lattice import g0k_from_ek
from triqs_tprf.lattice import gk_from_ek_sigma

from triqs_tprf.lattice import gr_from_gk
from triqs_tprf.lattice import grt_from_grw

from triqs_tprf.lattice import chi0_tr_from_grt_PH
from triqs_tprf.lattice import chi0_w0r_from_grt_PH
from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_w0r_from_chi_tr
from triqs_tprf.lattice import chi_wk_from_chi_wr

# ----------------------------------------------------------------------
def strip_sigma(nw, beta, sigma_in, debug=False):

    np.testing.assert_almost_equal(beta, sigma_in.mesh.beta)

    wmesh = MeshImFreq(beta, 'Fermion', n_max=nw)
    sigma = Gf(mesh=wmesh, target_shape=sigma_in.target_shape)

    for w in wmesh:
        index = w.linear_index + wmesh.first_index() # absolute index
        sigma[w] = sigma_in[Idx(index)]

    if debug:
        from pytriqs.plot.mpl_interface import oplot, plt
        oplot(p.Sigmalatt_iw)
        oplot(sigma,'x')
        plt.show(); exit()

    return sigma

# ----------------------------------------------------------------------
def bubble_setup(beta, mu, tb_lattice, nk, nw, sigma=None):

    print tprf_banner(), "\n"

    print 'beta  =', beta
    print 'mu    =', mu
    print 'sigma =', (not (sigma == None))

    norb = tb_lattice.NOrbitalsInUnitCell
    print 'nk    =', nk
    print 'nw    =', nw
    print 'norb  =', norb
    print

    ntau = 4 * nw
    ntot = np.prod(nk) * norb**4 + np.prod(nk) * (nw + ntau) * norb**2
    nbytes = ntot * np.complex128().nbytes
    ngb = nbytes / 1024.**3
    print 'Approx. Memory Utilization: %2.2f GB\n' % ngb
    
    periodization_matrix = np.diag(np.array(list(nk), dtype=np.int32))
    #print 'periodization_matrix =\n', periodization_matrix

    bz = BrillouinZone(tb_lattice.bl)
    bzmesh = MeshBrillouinZone(bz, periodization_matrix)

    print '--> ek'
    ek = ek_tb_dispersion_on_bzmesh(tb_lattice, bzmesh, bz)

    if sigma is None:
        print '--> g0k'
        wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)
        gwk = g0k_from_ek(mu=mu, ek=ek, mesh=wmesh)
    else:
        print '--> gk'
        sigma = strip_sigma(nw, beta, sigma)
        gwk = gk_from_ek_sigma(mu=mu, ek=ek, sigma=sigma)

    print '--> gr_from_gk (k->r)'
    gwr = gr_from_gk(gwk)
    del gwk

    print '--> grt_from_grw (w->tau)' 
    grt = grt_from_grw(gwr)
    del gwr

    if sigma is None:
        return grt
    else:
        return grt, sigma

# ----------------------------------------------------------------------
def chi0_w0k_tau_bubble(beta, mu, tb_lattice, nk, nw, sigma=None):

    if sigma is None:
        grt = bubble_setup(beta, mu, tb_lattice, nk, nw, sigma=sigma)
    else:
        grt, sigma_cut = bubble_setup(beta, mu, tb_lattice, nk, nw, sigma=sigma)
    
    print '--> chi0_w0r_from_grt_PH (bubble in tau & r)'
    chi0_wr = chi0_w0r_from_grt_PH(grt)
    del grt

    print '--> chi_wk_from_chi_wr (r->k)'
    chi0_wk = chi_wk_from_chi_wr(chi0_wr)
    del chi0_wr

    if sigma is None:
        return chi0_wk
    else:
        return chi0_wk, sigma_cut

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
def get_kidx_from_k_vec_relative(k_vec_rel, nk):
    kidx = np.array(np.round(k_vec_rel * np.array(nk)[None, :]), dtype=np.int)
    return kidx

# ----------------------------------------------------------------------
def get_k_components_from_k_vec(k_vec, nk):

    dim = 3

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
def get_abs_k_chi_interpolator(values, bzmesh, bz, extend_bz=[0]):

    k_mat = bz.units()
    k_vec = np.array([k.value for k in bzmesh])
    
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
def get_rel_k_chi_interpolator(values, bzmesh, bz, nk,
                               extend_boundary=True, interpolator='regular'):

    k_mat = bz.units()
    k_vec = np.array([k.value for k in bzmesh])

    k_vec_rel = get_relative_k_from_absolute(k_vec, bz.units())
    k_idx = get_kidx_from_k_vec_relative(k_vec_rel, nk)

    kx, ky, kz = get_k_components_from_k_vec(k_vec_rel, nk)

    if extend_boundary:
        values, k_vec_rel, (kx, ky, kz) = \
            extend_data_on_boundary(values, nk)
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
def extend_data_on_boundary(values, nk):

    # -- extended cube

    #nk = np.diag(periodization_matrix)
    nk = np.array(nk)

    # -- this adds points on the boundary

    # Add extra points in the positive index directions
    nk_ext = nk + 1
    coords = [ np.arange(0, n) for n in nk_ext ]

    # This adds two extra points in both positive and negative directions
    #nk_ext = nk + 4
    #coords = [ np.arange(-2, n+2) for n in nk ]
    
    Coords = np.meshgrid(*coords, indexing='ij')
    Coords_mod = [ np.mod(x, n) for x, n in zip(Coords, nk) ]

    values_ext = values.reshape(nk)[tuple([ X.flatten() for X in Coords_mod])]
    values_ext = values_ext.reshape(nk_ext)

    # -- compute kidx_ext

    k_idx_ext = np.array([ X.flatten() for X in Coords ]).T
    k_vec_rel_ext = np.array(k_idx_ext, dtype=np.float) / nk[None, :]
    kxe, kye, kze = get_k_components_from_k_vec(k_vec_rel_ext, nk_ext)

    return values_ext, k_vec_rel_ext, (kxe, kye, kze)
