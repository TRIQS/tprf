################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
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
import numpy as np

# ----------------------------------------------------------------------

import pytriqs.utility.mpi as mpi

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

from triqs_tprf.lattice import dynamical_screened_interaction_W_wk as cpp_dynamical_screened_interaction_W_wk
from triqs_tprf.lattice import \
    dynamical_screened_interaction_W_wk_from_generalized_susceptibility \
    as cpp_dynamical_screened_interaction_W_wk_from_generalized_susceptibility

from triqs_tprf.lattice import gw_sigma_wk_serial_fft as cpp_gw_sigma_wk_serial_fft
from triqs_tprf.lattice import gw_sigma_tr as cpp_gw_sigma_tr

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

    nw = len(g_wk.mesh.components[0]) / 2
    
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
def dynamical_screened_interaction_W_wk(PI_wk, V_k):
    return cpp_dynamical_screened_interaction_W_wk(PI_wk, V_k)

# ----------------------------------------------------------------------
def dynamical_screened_interaction_W_wk_from_generalized_susceptibility(chi_wk, V_k):
    return cpp_dynamical_screened_interaction_W_wk_from_generalized_susceptibility(chi_wk, V_k)

# ----------------------------------------------------------------------
def gw_sigma_wk(Wr_wk, g_wk, fft_flag=False):

    r""" GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator

    Fourier transforms the screened interaction and the single-particle
    Green's function to imagiary time and real space.

    .. math::
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

    .. math::
        W^{(r)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W^{(r)}_{abcd}(i\omega_n, \mathbf{k}) \right\}

    computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W^{(r)}_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    and transforms back to frequency and momentum

    .. math::
        \Sigma_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma_{ab}(\tau, \mathbf{r}) \right\}

    Parameters
    ----------

    V_k : TRIQS Green's function (rank 4) on a Brillouinzone mesh
          static bare interaction :math:`V_{abcd}(\mathbf{k})`

    Wr_wk : TRIQS Green's function (rank 4) on Matsubara and Brillouinzone meshes
            retarded screened interaction :math:`W^{(r)}_{abcd}(i\omega_n, \mathbf{k})`
    
    g_wk : TRIQS Green's function (rank 2) on Matsubara and Brillouinzone meshes
           single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

    Returns
    -------

    sigma_wk : TRIQS Green's function (rank 2) on Matsubara and Brillouinzone meshes
               GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`
    """

    if fft_flag:

        nw = len(g_wk.mesh.components[0]) / 2
        ntau = nw * 6 + 1
        
        mpi.report('g wk -> wr')
        g_wr = fourier_wk_to_wr(g_wk)
        mpi.report('g wr -> tr')
        g_tr = fourier_wr_to_tr(g_wr, nt=ntau)
        del g_wr

        mpi.report('W wk -> wr')
        Wr_wr = chi_wr_from_chi_wk(Wr_wk)
        mpi.report('W wr -> tr')
        Wr_tr = chi_tr_from_chi_wr(Wr_wr, ntau=ntau)
        del Wr_wr

        mpi.report('sigma tr')
        sigma_tr = cpp_gw_sigma_tr(Wr_tr, g_tr)
        del Wr_tr
        del g_tr

        mpi.report('sigma tr -> wr')
        sigma_wr = fourier_tr_to_wr(sigma_tr, nw=nw)

        del sigma_tr
        mpi.report('sigma wr -> wk')
        sigma_wk = fourier_wr_to_wk(sigma_wr)
        del sigma_wr

    else:
        sigma_wk = cpp_gw_sigma_wk_serial_fft(Wr_wk, g_wk)    
    
    return sigma_wk
