# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2018 by The Simons Foundation
# Copyright (C) 2020, S. Käser
# Copyright (C) 2023, Hugo U. R. Strand
# Authors: H. U.R. Strand, S. Käser 
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

import triqs.utility.mpi as mpi

from triqs.gf import Gf
from triqs.gf import Idx
from triqs.gf import MeshImFreq
from triqs.gf import MeshDLRImFreq
from triqs.gf import MeshProduct
from triqs.gf import MeshBrZone

from triqs.gf import make_gf_dlr

from triqs.lattice import BrillouinZone

# ----------------------------------------------------------------------

from triqs_tprf.logo import tprf_banner

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr

from triqs_tprf.lattice import chi0_tr_from_grt_PH
from triqs_tprf.lattice import chi0_wr_from_grt_PH
from triqs_tprf.lattice import chi0_w0r_from_grt_PH
from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_w0r_from_chi_tr
from triqs_tprf.lattice import chi_wk_from_chi_wr

from triqs_tprf.lattice import dlr_on_imfreq

# ----------------------------------------------------------------------
def add_fake_bosonic_mesh(gf, beta=None):
    """ Put a one value bosonic mesh as the first mesh argument of a 
    Green's function object.

    Parameters
    ----------
    gf : Gf,
         Green's function on some arbitrary mesh. If 'beta' is not given
         one mesh needs to be a 'MeshImFreq' to obtain a beta'
    beta : float, optional
           The inverse temperature used for the fake bosonic mesh.

    Returns
    -------
    gf_w : Gf,
           Green's function with an additional one value bosonic mesh
           on its first position.
    """
    mesh = gf.mesh
    if isinstance(mesh, MeshProduct):
        meshes = mesh.components
    else:
        meshes = (mesh,)

    # If beta is not given access it from a 'MeshImFreq' of the 'Gf'
    if not beta:
        betas = [mesh.beta for mesh in meshes if hasattr(mesh, "beta")]
        if len(betas) == 0:
            raise ValueError(
            "No 'beta' was given and the Green's function does not contain"
            " a 'MeshImFreq'")
        beta = betas[0]

    wmesh = MeshImFreq(beta, 'Boson', 1)
    mesh = (wmesh,) + meshes
    mesh = MeshProduct(*mesh)

    gf_w = Gf(mesh=mesh, target_shape=gf.target_shape)
    gf_w.data[0,...] = gf.data

    return gf_w

# ----------------------------------------------------------------------
def put_gf_on_mesh(g_in, wmesh):

    assert( len(wmesh) <= len(g_in.mesh) )
    
    g_out = Gf(mesh=wmesh, target_shape=g_in.target_shape)

    for w in wmesh:
        index = w.data_index + wmesh.first_index() # absolute index
        g_out[w] = g_in[Idx(index)]

    return g_out

# ----------------------------------------------------------------------
def strip_sigma(nw, beta, sigma_in, debug=False):

    np.testing.assert_almost_equal(beta, sigma_in.mesh.beta)

    wmesh = MeshImFreq(beta, 'Fermion', n_max=nw)
    sigma = Gf(mesh=wmesh, target_shape=sigma_in.target_shape)

    for w in wmesh:
        index = w.data_index + wmesh.first_index() # absolute index
        sigma[w] = sigma_in[Idx(index)]

    if debug:
        from triqs.plot.mpl_interface import oplot, plt
        oplot(p.Sigmalatt_iw)
        oplot(sigma,'x')
        plt.show(); exit()

    return sigma

# ----------------------------------------------------------------------
def bubble_setup(beta, mu, tb_lattice, nk, nw, sigma_w=None):

    print((tprf_banner(), "\n"))

    print(('beta  =', beta))
    print(('mu    =', mu))
    print(('sigma =', (not (sigma == None))))

    norb = tb_lattice.NOrbitalsInUnitCell
    print(('nk    =', nk))
    print(('nw    =', nw))
    print(('norb  =', norb))
    print()

    ntau = 4 * nw
    ntot = np.prod(nk) * norb**4 + np.prod(nk) * (nw + ntau) * norb**2
    nbytes = ntot * np.complex128().nbytes
    ngb = nbytes / 1024.**3
    print(('Approx. Memory Utilization: %2.2f GB\n' % ngb))
    
    bz = BrillouinZone(tb_lattice.bl)
    bzmesh = MeshBrZone(bz, nk)

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
def imtime_bubble_chi0_wk(g_wk, nw=1, save_memory=False, verbose=True):
    ncores = multiprocessing.cpu_count()

    wmesh, kmesh =  g_wk.mesh.components

    norb = g_wk.target_shape[0]
    beta = wmesh.beta
    nw_g = len(wmesh)
    nk = len(kmesh)

    ntau = 2 * nw_g

    is_dlr_mesh = type(wmesh) == MeshDLRImFreq
    if is_dlr_mesh:
        iw_values = np.fromiter(wmesh.values(), dtype=complex)
        is_symmetrized = np.allclose(iw_values, -iw_values[::-1])

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

    if verbose and mpi.is_master_node():
        print(tprf_banner(), "\n")
        print('beta  =', beta)
        print('nk    =', nk)
        print('nw    =', nw_g)
        print('norb  =', norb)
        print()
        print('Approx. Memory Utilization: %2.2f GB\n' % ngb)

    if verbose: mpi.report('--> fourier_wk_to_wr')
    g_wr = fourier_wk_to_wr(g_wk)
    del g_wk

    if verbose: mpi.report('--> fourier_wr_to_tr')
    g_tr = fourier_wr_to_tr(g_wr)
    del g_wr
    
    if nw == 1:
        if verbose: mpi.report('--> chi0_w0r_from_grt_PH (bubble in tau & r)')
        chi0_wr = chi0_w0r_from_grt_PH(g_tr)
        del g_tr
    else:
        if not save_memory:
            if verbose: mpi.report('--> chi0_tr_from_grt_PH (bubble in tau & r)')
            if is_dlr_mesh:
                chi0_tr = chi0_tr_from_grt_PH(g_tr, is_symmetrized)
            else:
                chi0_tr = chi0_tr_from_grt_PH(g_tr)
            del g_tr
            
            if verbose: mpi.report('--> chi_wr_from_chi_tr')
            chi0_wr = chi_wr_from_chi_tr(chi0_tr, nw=nw)
            del chi0_tr
        elif save_memory:
            chi0_wr = chi0_wr_from_grt_PH(g_tr, nw=nw)

        
    if verbose: mpi.report('--> chi_wk_from_chi_wr (r->k)')
    chi0_wk = chi_wk_from_chi_wr(chi0_wr)
    del chi0_wr

    return chi0_wk
# ----------------------------------------------------------------------
def imtime_bubble_chi0_ab_wk(ga_wk,gb_wk, nw=1, save_memory=False, verbose=True):
    ncores = multiprocessing.cpu_count()

    wmesh, kmesh =  g_wk.mesh.components

    norb = g_wk.target_shape[0]
    beta = wmesh.beta
    nw_g = len(wmesh)
    nk = len(kmesh)

    ntau = 2 * nw_g

    is_dlr_mesh = type(wmesh) == MeshDLRImFreq
    if is_dlr_mesh:
        iw_values = np.fromiter(wmesh.values(), dtype=complex)
        is_symmetrized = np.allclose(iw_values, -iw_values[::-1])

    
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

    if verbose and mpi.is_master_node():
        print(tprf_banner(), "\n")
        print('beta  =', beta)
        print('nk    =', nk)
        print('nw    =', nw_g)
        print('norb  =', norb)
        print()
        print('Approx. Memory Utilization: %2.2f GB\n' % ngb)

    if verbose: mpi.report('--> fourier_wk_to_wr')
    ga_wr = fourier_wk_to_wr(ga_wk)
    gb_wr = fourier_wk_to_wr(gb_wk)
    del ga_wk
    del gb_wk

    if verbose: mpi.report('--> fourier_wr_to_tr')
    ga_tr = fourier_wr_to_tr(ga_wr)
    gb_tr = fourier_wr_to_tr(gb_wr)
    del ga_wr
    del gb_wr
    
    if nw == 1:
        if verbose: mpi.report('--> chi0_w0r_from_grt_PH (bubble in tau & r)')
        chi0_wr = chi0_w0r_from_grt_ab_PH(ga_tr, gb_tr)
        del ga_tr
        del gb_tr
    else:
        if not save_memory:
            if verbose: mpi.report('--> chi0_tr_from_grt_PH (bubble in tau & r)')
            chi0_tr = chi0_tr_from_grt_ab_PH(ga_tr, gb_tr)
            del ga_tr
            del gb_tr
            
            if verbose: mpi.report('--> chi_wr_from_chi_tr')
            chi0_wr = chi_wr_from_chi_tr(chi0_tr, nw=nw)
            del chi0_tr
        elif save_memory:
            chi0_wr = chi0_wr_from_grt_PH(g_tr, nw=nw)

        
    if verbose: mpi.report('--> chi_wk_from_chi_wr (r->k)')
    chi0_wk = chi_wk_from_chi_wr(chi0_wr)
    del chi0_wr

    return chi0_wk
# ----------------------------------------------------------------------
def chi_contraction(chi, op1, op2):
    """Contract a susceptibility with two operators

    Parameters
    ----------
    chi : Gf,
          Susceptibility :math:`\chi(i\omega_n, \mathbf{k})`. The mesh attribute of
          the Gf must be a MeshProduct with the components (MeshImFreq, MeshBrZone)
          and its target_rank 4.
    op1, op2 : np.ndarray,
               Operators in matrix representation.

    Returns
    -------
    Gf,
    Susceptibility :math:`\chi(i\omega_n, \mathbf{k})`. With a target_rank of 0.
    """
    if chi.target_shape[:2] != op1.shape or chi.target_shape[2:] != op2.shape:
        raise ValueError('The shape of the operators %s and %s'%(op1.shape, op2.shape) +
                         ' must fit the shape of chi %s.'%(chi.target_shape,))

    chi_op1op2 = chi[0, 0, 0, 0].copy()
    chi_op1op2.data[:] = np.einsum('...abcd,ab,cd->...', chi.data, op1, op2)

    return chi_op1op2

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
    k_vec_rel = get_relative_k_from_absolute(k_vec, bz.units)
    
    ek.data[:] = tb_lattice.hopping(k_vec_rel.T).transpose(2, 0, 1)

    return ek

# ----------------------------------------------------------------------
def get_relative_k_from_absolute(k_vec, units):
    k_vec_rel = np.dot(np.linalg.inv(units).T, k_vec.T).T
    return k_vec_rel

# ----------------------------------------------------------------------
def get_kidx_from_k_vec_relative(k_vec_rel, nk):
    kidx = np.array(np.round(k_vec_rel * np.array(nk)[None, :]), dtype=np.int64)
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

    k_mat = bz.units
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

    k_mat = bz.units
    k_vec = np.array([k.value for k in bzmesh])

    k_vec_rel = get_relative_k_from_absolute(k_vec, bz.units)
    k_idx = get_kidx_from_k_vec_relative(k_vec_rel, nk)

    kx, ky, kz = get_k_components_from_k_vec(k_vec_rel, nk)

    if extend_boundary:
        values, k_vec_rel, (kx, ky, kz) = \
            extend_data_on_boundary(values, nk)
    else:
        values = values.reshape(nk)

    # -- select interpolator type
        
    if interpolator == 'regular':
        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator(
            (kx, ky, kz), values, fill_value=float('nan'), bounds_error=False)
    elif interpolator == 'nearest':
        from scipy.interpolate import NearestNDInterpolator        
        interp = NearestNDInterpolator(k_vec_rel, values.flatten())
    elif interpolator == 'linear':
        from scipy.interpolate import LinearNDInterpolator
        interp = LinearNDInterpolator(k_vec_rel, values.flatten(), fill_value=float('nan'))
    else:
        raise NotImplementedError
        
    return interp

# ----------------------------------------------------------------------
def extend_data_on_boundary(values, nk):

    # -- extended cube

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
def k_space_path(paths, num=100, bz=None, relative_coordinates=False):

    """ Generate an array of k-vectors along a path defined by a list of pairs of k-vectors

    Parameters
    ----------
    paths : list of pairs of three-vectors of floats
       List of pairs of k-vectors in reciprocal units to create a path in-between.
    num : int, default=100
       Number of k-vectors along each segment of the overall path
    bz : brillouin_zone, optional
       When a Brillouin Zone is passed, calculate distance in absolute units
    relative_coordinates : bool, optional
        Return k-vectors in reciprocal units. (Default `False`)

    Returns
    -------
    kvecs: numpy.ndarray [shape=(len(paths)*num,3)]
        Two-dimensional numpy array containing the path vectors (in reciprocal units) as rows
    dist: numpy.ndarray  [shape=(kvecs.shape[0])]
        One-dimensional numpy array containing, for each element in kvecs,
        the distance travelled along the path. Useful for plotting.
        If bz is provided, calculate the distance in absolute units.
        The distances for the relevant k-points in paths can be obtained with dist[num::num].
    ticks : numpy.ndarray [shape=(len(paths)+1)]
        Array with tick points, i.e. distances at the interfaces of path segments.
        Only returned if `return_ticks` is `True`.
    """
    
    # Put this back once TPRF and TRIQS agree on default output units
    #print("WARNING: triqs_tprf.lattice_utils.k_space_path has moved to triqs.lattice.utils.k_space_path")

    from triqs.lattice.utils import k_space_path
    return k_space_path(
        paths, num=num, bz=bz, relative_coordinates=relative_coordinates)


# ----------------------------------------------------------------------
def gf_tensor_to_matrix(g):
    
    N = g.target_shape[0] * g.target_shape[1] 
    M = g.target_shape[2] * g.target_shape[3] 
    g_mat = Gf(mesh=g.mesh, target_shape=[N, M])
    
    g_mat.data[:] = np.reshape(
        g.data, (g.data.shape[0], N, M))
    return g_mat


# ----------------------------------------------------------------------
def gf_matrix_to_tensor(g_mat, target_shape):
    
    g = Gf(mesh=g_mat.mesh, target_shape=target_shape)
    shape = [g_mat.data.shape[0]] + list(target_shape)
    g.data[:] = np.reshape(g_mat.data, shape)
    return g


# ----------------------------------------------------------------------
def pade_analytical_continuation_wk(
    g_wk, fmesh, n_points=32, freq_offset=0.05):

    """ Perform Pade analytical continuation of a lattice Green's function

    Parameters
    ----------

    g_wk : triqs.gf.Gf with mesh (MeshImFreq or MeshDLRImFreq, MeshBrZone)
        Lattice Green's function in Matsubara frequency and momentum space :math:`G(i\omega_n, \mathbf{k})`.
    fmesh : MeshReFreq
        Real frequency mesh to perform the analytical continuation to.
    n_points : int
        Number of Matsubara frequencies to use in the Pade fit.
    freq_offset : float
        Distance from the real axis used in the fit.

    Returns
    -------

    g_fk : triqs.gf.Gf with mesh (MeshReFreq, MeshBrZone)
        Analytically continued real-frequency lattice Green's function
    
    """
    
    wmesh = g_wk.mesh[0]
    kmesh = g_wk.mesh[1]

    g_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=g_wk.target_shape)

    for k in kmesh:
        g_f = g_fk[:, k]

        # Make_gf_dlr is broken on gf views (issue #913),
        # so make a copy of the data instead
        if type(wmesh) == MeshDLRImFreq:
            g_w = Gf(mesh=wmesh, target_shape=g_wk.target_shape)
            g_w.data[:] = g_wk[:, k].data[:]
        else:
            g_w = g_wk[:, k]

        if len(g_wk.target_shape) == 4:
            g_w = gf_tensor_to_matrix(g_w)
            g_f = gf_tensor_to_matrix(g_f)

        if type(g_w.mesh) == MeshDLRImFreq:
            g_c = make_gf_dlr(g_w)
            small_mesh = MeshImFreq(g_w.mesh.beta, g_w.mesh.statistic, n_points)
            g_w = dlr_on_imfreq(g_c, small_mesh)
            
        g_f.set_from_pade(g_w, n_points=n_points, freq_offset=freq_offset)

        if len(g_wk.target_shape) == 4:
            g_f = gf_matrix_to_tensor(g_f, g_wk.target_shape)
            g_fk[:, k] = g_f
        
    return g_fk
