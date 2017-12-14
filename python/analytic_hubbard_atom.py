################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2017 by Hugo U.R. Strand
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import *

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection

from triqs_tprf.hubbard_atom import chi_ph_magnetic
from triqs_tprf.hubbard_atom import gamma_ph_magnetic
from triqs_tprf.hubbard_atom import single_particle_greens_function

from triqs_tprf.linalg import inverse_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued

# ----------------------------------------------------------------------
def analytic_hubbard_atom(beta, U, nw, nwf, nwf_gf):

    d = ParameterCollection()
    d.beta, d.U, d.nw, d.nwf, d.nwf_gf = beta, U, nw, nwf, nwf_gf
    
    g_iw = single_particle_greens_function(beta=beta, U=U, nw=nwf_gf)
    d.G_iw = g_iw

    # make block gf of the single gf
    G_iw_block = BlockGf(name_list=['up', 'dn'], block_list=[g_iw, g_iw])
    g_mat = block_iw_AB_to_matrix_valued(G_iw_block)
    
    d.chi_m = chi_ph_magnetic(beta=beta, U=U, nw=nw, nwf=nwf)
    d.chi0_m = chi0_from_gg2_PH(g_mat, d.chi_m)

    # -- Numeric vertex from BSE
    d.gamma_m_num = inverse_PH(d.chi0_m) - inverse_PH(d.chi_m)

    # -- Analytic vertex
    d.gamma_m = gamma_ph_magnetic(beta=beta, U=U, nw=nw, nwf=nwf)

    # -- Analytic magnetization expecation value
    # -- and static susceptibility

    d.Z = 2. + 2*np.exp(-beta * 0.5 * U)
    d.m2 = 0.25 * (2 / d.Z)
    d.chi_m_static = 2. * beta * d.m2
    
    d.label = r'Analytic'
    
    return d

# ----------------------------------------------------------------------
