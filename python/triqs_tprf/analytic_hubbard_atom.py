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

    r""" Compute dynamical response functions for the Hubbard atom at half filling.

    This function returns an object that contains the single-particle
    Greens function :math:`G(\omega)`, the magnetic two-particle generalized susceptibility
    :math:`\chi_m(\omega, \nu, \nu')`, and the corresponding bare bubble 
    :math:`\chi^{(0)}_m(\omega, \nu, \nu')`, and the magnetic vertex function
    :math:`\Gamma_m(\omega, \nu, \nu')`.

    This is implemented using analytical formulas from 
    Thunstrom et al. [PRB 98, 235107 (2018)]
    please cite the paper if you use this function!

    In particular this is one exact solution to the Bethe-Salpeter
    equation, that is the infinite matrix inverse problem:

    .. math::
        \Gamma_m = [\chi^{(0)}_m]^{-1} - \chi_m^{-1}

    Parameters
    ----------

    beta : float
        Inverse temperature.

    U : float
        Hubbard U interaction parameter.

    nw : int
        Number of bosonic Matsubara frequencies 
        in the computed two-particle response functions.

    nwf : int
        Number of fermionic Matsubara frequencies
        in the computed two-particle response functions.

    nwf_gf : int
        Number of fermionic Matsubara frequencies
        in the computed single-particle Greens function.

    Returns
    -------

    p : ParameterCollection
        Object containing all the response functions and some other
        observables, `p.G_iw`, `p.chi_m`, `p.chi0_m`, 
        `p.gamma_m`, `p.Z`, `p.m2`, `p.chi_m_static`.

    """
    
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
