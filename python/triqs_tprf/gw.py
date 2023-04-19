################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
# Copyright (C) 2023 by Hugo U.R. Strand
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

import itertools
import numpy as np

# ----------------------------------------------------------------------

import triqs.utility.mpi as mpi

# ----------------------------------------------------------------------

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_wk
from triqs_tprf.lattice import fourier_wr_to_tr
from triqs_tprf.lattice import fourier_tr_to_wr

from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_tr_from_chi_wr
from triqs_tprf.lattice import chi_wk_from_chi_wr
from triqs_tprf.lattice import chi_wr_from_chi_wk

from triqs_tprf.lattice import chi0_tr_from_grt_PH

from triqs_tprf.lattice import dynamical_screened_interaction_W
from triqs_tprf.lattice import \
    dynamical_screened_interaction_W_from_generalized_susceptibility

from triqs_tprf.lattice import gw_sigma, g0w_sigma
from triqs_tprf.lattice import lindhard_chi00

from triqs_tprf.OperatorUtils import quartic_tensor_from_operator

# ----------------------------------------------------------------------
def bubble_PI_wk(g_wk):

    r""" Compute the particle-hole bubble from the single particle lattice Green's function

    .. math::
        \Pi_{abcd}(i\omega_n, k) = - 
        \mathcal{F}_{\tau, \mathbf{r} \rightarrow i\omega_n, \mathbf{k}}
        \left\{ G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r}) \right\}
    
    Parameters
    ----------

    g_wk : TRIQS Green's function (rank 2) on Matsubara and Brillouinzone meshes
           Single particle lattice Green's function.

    Returns
    -------

    PI_wk : TRIQS Green's function (rank 4) on Matsubara and Brillouinzone meshes
            Particle hole bubble.
    """

    nw = len(g_wk.mesh.components[0]) // 2
    
    g_wr = fourier_wk_to_wr(g_wk)
    g_tr = fourier_wr_to_tr(g_wr)
    del g_wr
    PI_tr = chi0_tr_from_grt_PH(g_tr)
    del g_tr
    PI_wr = chi_wr_from_chi_tr(PI_tr, nw=nw)
    del PI_tr
    PI_wk = chi_wk_from_chi_wr(PI_wr)
    del PI_wr

    return PI_wk


# ----------------------------------------------------------------------
def get_gw_tensor(H_int, fundamental_operators):

    """ Takes a TRIQS operator object and extracts the quartic terms
    and returns the corresponding antisymmetrized quartic tensor in
    vertex index order, i.e., cc+cc+. """
    
    U_abcd = quartic_tensor_from_operator(H_int, fundamental_operators)
    
    # -- Group in c^+cc^+c ( from c^+c^+cc )
    U_abcd = np.ascontiguousarray(np.transpose(U_abcd, (0, 3, 1, 2)))    

    # -- Density-density permutations
    U_abcd = (U_abcd + np.transpose(U_abcd, (2,3,0,1)))

    return U_abcd
