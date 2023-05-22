# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 H. U.R. Strand
# Authors: H. U.R. Strand
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

from triqs.gf import Gf, MeshProduct, Idx

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi0r_from_gr_PH
from triqs_tprf.lattice import chi0q_from_chi0r
from triqs_tprf.lattice import chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH

from triqs_tprf.linalg import product_PH, inverse_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH


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

    chi_wnn = chi_from_gg2_PH(g_w, g2_wnn)
    chi0_wnn = chi0_from_gg2_PH(g_w, g2_wnn)

    g2_conn_wnn = chi_wnn - chi0_wnn

    inv_chi0_wnn = inverse_PH(chi0_wnn)
    F_wnn = product_PH(inv_chi0_wnn, product_PH(g2_conn_wnn, inv_chi0_wnn))
    
    return F_wnn


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

    bmesh = F_wnn.mesh[0]
    fmesh = F_wnn.mesh[1]

    assert( len(fmesh) <= len(g_wk.mesh[0]) )
    assert( len(bmesh) <= len(L_wn.mesh[0]) )
    assert( len(fmesh) <= len(L_wn.mesh[1]) )
    assert( len(bmesh) <= len(chi_imp_w.mesh) )
    
    nw = (len(bmesh) + 1) // 2
    nn = len(fmesh) // 2
    
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
