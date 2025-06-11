# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 H. U.R. Strand and E G.C.P. van Loon
# Authors: H. U.R. Strand and E. G.C.P. van Loon
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

import numpy as np

from triqs.gf import Gf, MeshProduct, Idx

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi0r_from_gr_PH
from triqs_tprf.lattice import chi0q_from_chi0r
from triqs_tprf.lattice import chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH
from triqs_tprf.lattice import chiq_sum_nu_from_e_k_sigma_w_F_nn_and_L_n_PH
from triqs_tprf.lattice import lattice_dyson_g_w

from triqs_tprf.linalg import product_PH, inverse_PH
from triqs_tprf.linalg_utils import gf_matrix_from_tensor, gf_tensor_from_matrix
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH

from triqs_tprf.bse import get_chi0_nk_at_specific_w
from triqs_tprf.lattice_utils import add_fake_bosonic_mesh


def impurity_reducible_vertex_F(g_w, g2_wnn):

    r"""Compute the impurity reducible vertex function 
    :math:`F_{abcd}(\omega, \nu, \nu')`.

    Computes:

    .. math::
       F_{abcd}(\omega, \nu, \nu') =  [\chi^{(0)}]^{-1} (\chi - \chi^{(0)} ) [\chi^{(0)}]^{-1} 

    where the inverses are taken in the particle-hole channel pairing
    of fermionic frequencies :math:`\nu` and :math:`\nu'` and orbital
    indices.

    Parameters
    ----------

    g_w : Single particle Green's function
          :math:`G_{ab}(\nu)`
    g2_wnn : Two-particle Green's function
             :math:`G^{(2)}_{abcd}(\omega, \nu, \nu')`

    Returns
    -------

    F_wnn : Particle-hole reducible vertex function 
            :math:`F_{abcd}(\omega, \nu, \nu')`
    """

    assert( len(g_w.target_shape) == 2 )
    assert( len(g2_wnn.target_shape) == 4 )

    norb = g_w.target_shape[0]

    assert( (np.array(g_w.target_shape) == norb).all() )
    assert( (np.array(g2_wnn.target_shape) == norb).all() )

    chi_wnn = chi_from_gg2_PH(g_w, g2_wnn)
    chi0_wnn = chi0_from_gg2_PH(g_w, g2_wnn)

    g2_conn_wnn = chi_wnn - chi0_wnn

    inv_chi0_wnn = inverse_PH(chi0_wnn)
    F_wnn = product_PH(inv_chi0_wnn, product_PH(g2_conn_wnn, inv_chi0_wnn))
    
    return F_wnn


def impurity_reducible_vertex_F_nn(w, g_w, g2_nn):

    r"""Compute the impurity reducible vertex function 
    :math:`F_{abcd}(\omega, \nu, \nu')` at fixed bosonic frequency :math:`\omega`.
    
    Computes:

    .. math::
       F_{abcd}(\omega, \nu, \nu') =  [\chi^{(0)}]^{-1} (\chi - \chi^{(0)} ) [\chi^{(0)}]^{-1} 

    where the inverses are taken in the particle-hole channel pairing
    of fermionic frequencies :math:`\nu` and :math:`\nu'` and orbital
    indices.

    Parameters
    ----------

    w : Bosonic Matsubara Frequency
    g_w : Single particle Green's function
          :math:`G_{ab}(\nu)`
    g2_nn : Two-particle Green's function
            :math:`G^{(2)}_{abcd}(\nu, \nu')`

    Returns
    -------

    F_nn : Particle-hole reducible vertex function 
           :math:`F_{abcd}(\omega, \nu, \nu')` at the given bosonic frequency :math:`\omega`.
    """
    
    fmesh = g2_nn.mesh[0]
    beta = fmesh.beta
    
    chi0_n = Gf(mesh=fmesh, target_shape=g2_nn.target_shape)    
    for n in fmesh:
        chi0_n[n] = -beta * np.einsum('da,bc->abcd', g_w(n).data, g_w(w + n).data)

    g2c_nn = g2_nn.copy()

    for n in fmesh:
        g2c_nn[n, n] -= chi0_n[n]

    if w.index == 0:
        # Correction
        for n1, n2 in g2c_nn.mesh:
            g2c_nn[n1, n2] -= beta * np.einsum('ba,dc->abcd', g_w(n1).data, g_w(n2).data)

    chi0_n_mat = gf_matrix_from_tensor(chi0_n)
    g2c_nn_mat = gf_matrix_from_tensor(g2c_nn)

    chi0_n_mat_inv = chi0_n_mat.copy()
    chi0_n_mat_inv.data[:] = np.linalg.inv(chi0_n_mat.data)
        
    F_nn_mat = g2c_nn_mat.copy()
    for n1, n2 in F_nn_mat.mesh:
        F_nn_mat[n1, n2] = np.matmul(
            chi0_n_mat_inv[n1].data, np.matmul(g2c_nn_mat[n1, n2].data, chi0_n_mat_inv[n2].data))

    F_nn = gf_tensor_from_matrix(F_nn_mat)
        
    return F_nn


def solve_lattice_dbse(g_wk, F_wnn, L_wn, chi_imp_w):

    r""" Compute the generalized lattice susceptibility 
    :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, \omega_n)` using the dual Bethe-Salpeter 
    equation (DBSE).

    Parameters
    ----------

    g_wk : Gf,
           Single-particle Green's function :math:`G_{a\bar{b}}(i\nu_n, \mathbf{k})`.
    F_wnn : Gf,
                Local particle-hole reducible vertex function 
                :math:`F_{a\bar{b}c\bar{d}}(i\omega_n, i\nu_n, i\nu_n')`.
    L_wn : Gf,
                Local particle-hole reducible triangle vertex function 
                :math:`L_{a\bar{b}c\bar{d}}(i\omega_n, i\nu_n)`.
    chi_imp_w : Gf,
                Generalized DMFT impurity susceptibility
                :math:`\chi_{a\bar{b}c\bar{d}}(i\omega_n)`.

    Returns
    -------
    chi_kw : Gf,
             Generalized lattice susceptibility 
             :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`.
    """

    # -- Check mesh sizes

    bmesh = F_wnn.mesh[0]
    fmesh = F_wnn.mesh[1]

    assert( len(fmesh) <= len(g_wk.mesh[0]) )
    assert( len(bmesh) <= len(L_wn.mesh[0]) )
    assert( len(fmesh) <= len(L_wn.mesh[1]) )
    assert( len(bmesh) <= len(chi_imp_w.mesh) )
    
    nw = (len(bmesh) + 1) // 2
    nn = len(fmesh) // 2

    # -- Check target_shape(s)

    assert( len(g_wk.target_shape) == 2 )
    assert( len(F_wnn.target_shape) == 4 )
    assert( len(L_wn.target_shape) == 4 )
    assert( len(chi_imp_w.target_shape) == 4 )

    norb = g_wk.target_shape[0]

    assert( (np.array(g_wk.target_shape) == norb).all() )
    assert( (np.array(F_wnn.target_shape) == norb).all() )
    assert( (np.array(L_wn.target_shape) == norb).all() )
    assert( (np.array(chi_imp_w.target_shape) == norb).all() )

    
    print('--> g_nonlocal_wr')
    # -- Remove local gf component (at r = 0)
    g_nonlocal_wr = fourier_wk_to_wr(g_wk)
    g_nonlocal_wr[:, Idx(0, 0, 0)] = 0.

    print('--> chi0_nonlocal_wnr')
    chi0_nonlocal_wnr = chi0r_from_gr_PH(nw=nw, nn=nn, g_nr=g_nonlocal_wr)

    del g_nonlocal_wr
    
    print('--> chi0_nonlocal_wnk')
    chi0_nonlocal_wnk = chi0q_from_chi0r(chi0_nonlocal_wnr)

    del chi0_nonlocal_wnr
    
    print('--> Resize L_wn')
    L_resize_wn = Gf(indices=L_wn.indices, mesh=MeshProduct(bmesh, fmesh))
    for w, n in L_resize_wn.mesh:
        L_resize_wn[w, n] = L_wn[Idx(w.index), Idx(n.index)]    
        
    print('--> DBSE chi_kw')
    chi_kw = chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH(
        chi0_nonlocal_wnk, F_wnn, L_resize_wn)
    
    for w in bmesh:
        chi_kw[:, w].data[:] += chi_imp_w[Idx(w.index)].data

    del chi0_nonlocal_wnk
    del L_resize_wn
    
    return chi_kw


def solve_lattice_dbse_lomem_w(g_wk, F_wnn, L_wn, chi_imp_w):

    r""" Compute the generalized lattice susceptibility 
    :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, \omega_n)` using the dual Bethe-Salpeter 
    equation (DBSE).

    This implementation loops over individual bosonic frequencies, requiring less memory
    but slightly more computations.

    Parameters
    ----------

    g_wk : Gf,
           Single-particle Green's function :math:`G_{a\bar{b}}(i\nu_n, \mathbf{k})`.
    F_wnn : Gf,
                Local particle-hole reducible vertex function 
                :math:`F_{a\bar{b}c\bar{d}}(i\omega_n, i\nu_n, i\nu_n')`.
    L_wn : Gf,
                Local particle-hole reducible triangle vertex function 
                :math:`L_{a\bar{b}c\bar{d}}(i\omega_n, i\nu_n)`.
    chi_imp_w : Gf,
                Generalized DMFT impurity susceptibility
                :math:`\chi_{a\bar{b}c\bar{d}}(i\omega_n)`.

    Returns
    -------
    chi_kw : Gf,
             Generalized lattice susceptibility 
             :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`.
    """

    # -- Check mesh sizes

    bmesh = F_wnn.mesh[0]
    fmesh = F_wnn.mesh[1]

    assert( len(fmesh) <= len(g_wk.mesh[0]) )
    assert( len(bmesh) <= len(L_wn.mesh[0]) )
    assert( len(fmesh) <= len(L_wn.mesh[1]) )
    assert( len(bmesh) <= len(chi_imp_w.mesh) )
    
    nw = (len(bmesh) + 1) // 2
    nn = len(fmesh) // 2

    # -- Check target_shape(s)

    assert( len(g_wk.target_shape) == 2 )
    assert( len(F_wnn.target_shape) == 4 )
    assert( len(L_wn.target_shape) == 4 )
    assert( len(chi_imp_w.target_shape) == 4 )

    norb = g_wk.target_shape[0]

    assert( (np.array(g_wk.target_shape) == norb).all() )
    assert( (np.array(F_wnn.target_shape) == norb).all() )
    assert( (np.array(L_wn.target_shape) == norb).all() )
    assert( (np.array(chi_imp_w.target_shape) == norb).all() )

    
    L_resize_wn = Gf(mesh=MeshProduct(bmesh, fmesh), indices=L_wn.indices)
    for w, n in L_resize_wn.mesh:
        L_resize_wn[w, n] = L_wn[Idx(w.index), Idx(n.index)]

    kmesh = g_wk.mesh[1]
    chi_kw = Gf(mesh=MeshProduct(kmesh, bmesh), indices=F_wnn.indices)
    
    for W in bmesh:
        
        print('-'*72)
        print(f'DBSE: Low memory calc at bosonic frequency index {W.index}.')
        print(f'DBSE: {W}')
        print('-'*72)

        idx = Idx(W.index)
        F_Wnn = add_fake_bosonic_mesh(F_wnn[idx, :, :])
        L_resize_Wn = add_fake_bosonic_mesh(L_resize_wn[idx, :])

        print('--> chi0_nonlocal_Wnk')
        chi0_nonlocal_Wnk = add_fake_bosonic_mesh(get_chi0_nk_at_specific_w(
            g_wk, nw_index=W.index, nwf=nn, g_nonlocal=True))

        print('--> DBSE chi_kW')
        chi_kW = chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH(
            chi0_nonlocal_Wnk, F_Wnn, L_resize_Wn)

        chi_kw[:, idx] = chi_kW[:, Idx(0)] + chi_imp_w[idx].data
            
        del chi0_nonlocal_Wnk
        del chi_kW

    del L_resize_wn

    return chi_kw


def solve_lattice_dbse_lomem_kw(mu, e_k, sigma_w, F_wnn, L_wn, chi_imp_w):

    r""" Compute the generalized lattice susceptibility 
    :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, \omega_n)` using the dual Bethe-Salpeter 
    equation (DBSE).

    This implementation loops over individual bosonic frequencies and mometa,
    requiring less memory but slightly more computations.

    Parameters
    ----------

    mu : double,
                Chemical potential
    e_k : Gf,
                Lattice dispersion :math:`\epsilon_{a\bar{b}}(\mathbf{k})`.
    sigma_w : Gf,
                Local self-energy :math:`\Sigma_{a\bar{b}}(i\nu_n)`.    
    F_wnn : Gf,
                Local particle-hole reducible vertex function 
                :math:`F_{a\bar{b}c\bar{d}}(i\omega_n, i\nu_n, i\nu_n')`.
    L_wn : Gf,
                Local particle-hole reducible triangle vertex function 
                :math:`L_{a\bar{b}c\bar{d}}(i\omega_n, i\nu_n)`.
    chi_imp_w : Gf,
                Generalized DMFT impurity susceptibility
                :math:`\chi_{a\bar{b}c\bar{d}}(i\omega_n)`.

    Returns
    -------
    chi_kw : Gf,
             Generalized lattice susceptibility 
             :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`.
    """

    # -- Check mesh sizes

    bmesh = F_wnn.mesh[0]
    fmesh = F_wnn.mesh[1]

    assert( len(fmesh) <= len(sigma_w.mesh) )
    assert( len(bmesh) <= len(L_wn.mesh[0]) )
    assert( len(fmesh) <= len(L_wn.mesh[1]) )
    assert( len(bmesh) <= len(chi_imp_w.mesh) )
    
    nw = (len(bmesh) + 1) // 2
    nn = len(fmesh) // 2

    # -- Check target_shape(s)

    assert( len(sigma_w.target_shape) == 2 )
    assert( len(F_wnn.target_shape) == 4 )
    assert( len(L_wn.target_shape) == 4 )
    assert( len(chi_imp_w.target_shape) == 4 )

    norb = sigma_w.target_shape[0]

    assert( (np.array(sigma_w.target_shape) == norb).all() )
    assert( (np.array(F_wnn.target_shape) == norb).all() )
    assert( (np.array(L_wn.target_shape) == norb).all() )
    assert( (np.array(chi_imp_w.target_shape) == norb).all() )
    
    
    L_resize_wn = Gf(mesh=MeshProduct(bmesh, fmesh), indices=L_wn.indices)
    for w, n in L_resize_wn.mesh:
        L_resize_wn[w, n] = L_wn[Idx(w.index), Idx(n.index)]

    g_loc_w = lattice_dyson_g_w(mu, e_k, sigma_w)

    kmesh = e_k.mesh
    chi_kw = Gf(mesh=MeshProduct(kmesh, bmesh), indices=F_wnn.indices)

    mem = np.prod(chi_kw.data.shape) * 128 / 8
    print(f'chi_kw.data.shape = {chi_kw.data.shape}')
    print(f'Memory estimate: {mem / 1024**3} GB')
    
    chi_kw.data[:] = 0 # Trigger full alloc
    
    import itertools
    for W, Q in itertools.product(bmesh, kmesh):
        
        print('-'*72)
        print(f'DBSE: Low memory calc at bosonic frequency index {W.index} and momentum index {Q.index}.')
        print(f'DBSE: {W}')
        print(f'DBSE: {Q}')
        print('-'*72)

        idx = Idx(W.index)

        F_nn = F_wnn[W, :, :]
        L_resize_n = L_resize_wn[W, :]
        
        chi_QW = chiq_sum_nu_from_e_k_sigma_w_F_nn_and_L_n_PH(
            mu, e_k, sigma_w, g_loc_w, F_nn, L_resize_n, W.data_index, Q.data_index, bmesh)

        chi_kw[Q, idx] = chi_QW + chi_imp_w[idx].data
            
        del chi_QW

    del L_resize_wn

    return chi_kw
