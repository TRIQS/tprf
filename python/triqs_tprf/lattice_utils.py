
################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2018 by The Simons Foundation
# Author: H. U.R. Strand
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import itertools
import multiprocessing
import numpy as np

# ----------------------------------------------------------------------

import pytriqs.utility.mpi as mpi

from pytriqs.gf import Gf
from pytriqs.gf import Idx
from pytriqs.gf import MeshImFreq
from pytriqs.gf import MeshProduct
from pytriqs.gf import MeshBrillouinZone

from pytriqs.lattice import BrillouinZone

# ----------------------------------------------------------------------

from triqs_tprf.logo import tprf_banner

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr

from triqs_tprf.lattice import chi0_tr_from_grt_PH
from triqs_tprf.lattice import chi0_w0r_from_grt_PH
from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_w0r_from_chi_tr
from triqs_tprf.lattice import chi_wk_from_chi_wr

# ----------------------------------------------------------------------
def put_gf_on_mesh(g_in, wmesh):

    assert( len(wmesh) <= len(g_in.mesh) )
    
    g_out = Gf(mesh=wmesh, target_shape=g_in.target_shape)

    for w in wmesh:
        index = w.linear_index + wmesh.first_index() # absolute index
        g_out[w] = g_in[Idx(index)]

    return g_out

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
def bubble_setup(beta, mu, tb_lattice, nk, nw, sigma_w=None):

    print(tprf_banner(), "\n")

    print('beta  =', beta)
    print('mu    =', mu)
    print('sigma =', (not (sigma == None)))

    norb = tb_lattice.NOrbitalsInUnitCell
    print('nk    =', nk)
    print('nw    =', nw)
    print('norb  =', norb)
    print()

    ntau = 4 * nw
    ntot = np.prod(nk) * norb**4 + np.prod(nk) * (nw + ntau) * norb**2
    nbytes = ntot * np.complex128().nbytes
    ngb = nbytes / 1024.**3
    print('Approx. Memory Utilization: %2.2f GB\n' % ngb)
    
    periodization_matrix = np.diag(np.array(list(nk), dtype=np.int32))
    #print 'periodization_matrix =\n', periodization_matrix

    bz = BrillouinZone(tb_lattice.bl)
    bzmesh = MeshBrillouinZone(bz, periodization_matrix)

    print('--> ek')
    e_k = ek_tb_dispersion_on_bzmesh(tb_lattice, bzmesh, bz)

    if sigma is None:
        print('--> g0k')
        wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)
        g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
    else:
        print('--> gk')
        sigma_w = strip_sigma(nw, beta, sigma)
        g_wk = lattice_dyson_g_wk(mu=mu, e_k=e_k, sigma_w=sigma_w)

    print('--> gr_from_gk (k->r)')
    g_wr = fourier_wk_to_wr(g_wk)
    del g_wk

    print('--> grt_from_grw (w->tau)') 
    g_tr = fourier_wr_to_tr(g_wr)
    del g_wr

    if sigma is None:
        return g_tr
    else:
        return g_tr, sigma_w

# ----------------------------------------------------------------------
def imtime_bubble_chi0_wk(g_wk, nw=1):
    ncores = multiprocessing.cpu_count()

    wmesh, kmesh =  g_wk.mesh.components

    norb = g_wk.target_shape[0]
    beta = wmesh.beta
    nw_g = len(wmesh)
    nk = len(kmesh)

    ntau = 2 * nw_g

    # -- Memory Approximation

    ng_tr = ntau * np.prod(nk) * norb**2 # storing G(tau, r)
    ng_wr = nw_g * np.prod(nk) * norb**2 # storing G(w, r)
    ng_t = ntau * norb**2 # storing G(tau)

    nchi_tr = ntau * np.prod(nk) * norb**4 # storing \chi(tau, r)
    nchi_wr = nw * np.prod(nk) * norb**4 # storing \chi(w, r)
    nchi_t = ntau * norb**4 # storing \chi(tau)
    nchi_w = nw * norb**4 # storing \chi(w)
    nchi_r = np.prod(nk) * norb**4 # storing \chi(r)

    if nw == 1:
        ntot_case_1 = ng_tr + ng_wr
        ntot_case_2 = ng_tr + nchi_wr + ncores*(nchi_t + 2*ng_t)
        ntot_case_3 = 4 * nchi_wr

        ntot = max(ntot_case_1, ntot_case_2, ntot_case_3)

    else:
        ntot_case_1 = ng_tr + nchi_tr + ncores*(nchi_t + 2*ng_t)
        ntot_case_2 = nchi_tr + nchi_wr + ncores*(nchi_w + nchi_t)
    
        ntot = max(ntot_case_1, ntot_case_2)

    nbytes = ntot * np.complex128().nbytes
    ngb = nbytes / 1024.**3

    if mpi.is_master_node():
        print(tprf_banner(), "\n")
        print('beta  =', beta)
        print('nk    =', nk)
        print('nw    =', nw_g)
        print('norb  =', norb)
        print()
        print('Approx. Memory Utilization: %2.2f GB\n' % ngb)

    mpi.report('--> fourier_wk_to_wr')
    g_wr = fourier_wk_to_wr(g_wk)
    del g_wk

    mpi.report('--> fourier_wr_to_tr')
    g_tr = fourier_wr_to_tr(g_wr)
    del g_wr
    
    if nw == 1:
        mpi.report('--> chi0_w0r_from_grt_PH (bubble in tau & r)')
        chi0_wr = chi0_w0r_from_grt_PH(g_tr)
        del g_tr
    else:
        mpi.report('--> chi0_tr_from_grt_PH (bubble in tau & r)')
        chi0_tr = chi0_tr_from_grt_PH(g_tr)
        del g_tr
        
        mpi.report('--> chi_wr_from_chi_tr')
        chi0_wr = chi_wr_from_chi_tr(chi0_tr, nw=nw)
        del chi0_tr
        
    mpi.report('--> chi_wk_from_chi_wr (r->k)')
    chi0_wk = chi_wk_from_chi_wr(chi0_wr)
    del chi0_wr

    return chi0_wk

# ----------------------------------------------------------------------
def chi0_w0k_tau_bubble(beta, mu, tb_lattice, nk, nw, sigma_w=None):

    if sigma_w is None:
        g_tr = bubble_setup(beta, mu, tb_lattice, nk, nw)
    else:
        g_tr, sigma_w_cut = bubble_setup(beta, mu, tb_lattice, nk, nw, sigma_w=sigma_w)
    
    print('--> chi0_w0r_from_grt_PH (bubble in tau & r)')
    chi0_wr = chi0_w0r_from_grt_PH(g_tr)
    del grt

    print('--> chi_wk_from_chi_wr (r->k)')
    chi0_wk = chi_wk_from_chi_wr(chi0_wr)
    del chi0_wr

    if sigma is None:
        return chi0_wk
    else:
        return chi0_wk, sigma_w_cut

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
    for axis in range(dim):
        cut = [0]*dim
        cut[axis] = slice(None)
        cut = [axis] + cut
        k_out.append(k_vec[tuple(cut)])
        
    return tuple(k_out)

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

    from scipy.interpolate import LinearNDInterpolator
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
        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator(
            (kx, ky, kz), values, fill_value=float('nan'), bounds_error=False)
    elif interpolator is 'nearest':
        from scipy.interpolate import NearestNDInterpolator        
        interp = NearestNDInterpolator(k_vec_rel, values.flatten())
    elif interpolator is 'linear':
        from scipy.interpolate import LinearNDInterpolator
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

# ----------------------------------------------------------------------
def k_space_path(paths, num=100, bz=None):

    """ High symmetry path k-vector generator.

    Input:
    paths : list of tuples of pairs of 3-vectors of k-points to
    make the path in between 
    num (optional) : number of k-vectors along each segment of the path
    bz (optional) : Brillouin zone, used to rescale from reative to absolute k-space lengths

    Returns:
    k_vecs: ndarray.shape = (n_k, 3) with all k-vectors. 
    k_plot: ndarray.shape = (n_k) one dimensional vector for plotting
    K_plot: ndarray.shape = (n_paths) positions of the start and end of each path 
    """

    if bz is None:
        cell = np.eye(3)
    else:
        cell = bz.units()
    
    k_vecs = []

    for path in paths:
        ki, kf = path
        x = np.linspace(0., 1., num=num)[:, None]
        k_vec = (1. - x) * ki[None, :] + x * kf[None, :]

        k_vecs.append(k_vec)

    def rel_to_abs(k_vec, cell):
        return np.einsum('ba,ib->ia', cell, k_vec)

    k_vec = k_vecs[0]
    k_vec_abs = rel_to_abs(k_vec, cell)
    k_plot = np.linalg.norm(k_vec_abs - k_vec_abs[0][None, :], axis=1)

    K_plot = [0.]
    for kidx, k_vec in enumerate(k_vecs[1:]):
        k_vec_abs = rel_to_abs(k_vec, cell)
        k_plot_new = np.linalg.norm(k_vec_abs - k_vec_abs[0][None, :], axis=1) + k_plot[-1]
        K_plot.append(k_plot[-1])
        k_plot = np.concatenate((k_plot, k_plot_new))

    K_plot.append(k_plot[-1])
    K_plot = np.array(K_plot)
    k_vecs = np.vstack(k_vecs)
    
    return k_vecs, k_plot, K_plot
