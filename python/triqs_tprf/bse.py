
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

import numpy as np

# ----------------------------------------------------------------------

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.gf import MeshImFreq, MeshProduct, Gf, Idx

# ----------------------------------------------------------------------

from triqs_tprf.logo import tprf_banner
from triqs_tprf.linalg import inverse_PH

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi0r_from_gr_PH
from triqs_tprf.lattice import chi0r_from_gr_PH_nompi
from triqs_tprf.lattice import chi0q_from_chi0r
from triqs_tprf.lattice import chi0q_sum_nu, chi0q_sum_nu_tail_corr_PH
from triqs_tprf.lattice import chiq_sum_nu_from_chi0q_and_gamma_PH

# ----------------------------------------------------------------------
def solve_local_bse(chi0_wnn, chi_wnn):

    r"""Solve the Bethe-Salpeter equation for the local vertex function 
    :math:`\Gamma_{abcd}(\omega, \nu, \nu')`.

    Computes:

    .. math::
       \Gamma_{abcd}(\omega, \nu, \nu') = [\chi^{(0)}]^{-1} - \chi^{-1}

    where the inverses are taken in the particle-hole channel pairing
    of fermionic frequencies :math:`\nu` and :math:`\nu'` and orbital
    indices.

    Parameters
    ----------

    chi0_wnn : Gerealized local bubble susceptibility 
               :math:`\chi^{(0)}_{abcd}(\omega, \nu, \nu')`
    chi_wnn : Generalized local susceptibility 
              :math:`\chi_{abcd}(\omega, \nu, \nu')`

    Returns
    -------

    gamma_wnn : Particle-hole vertex function 
                :math:`\Gamma_{abcd}(\omega, \nu, \nu')`
    """

    gamma_wnn = inverse_PH(chi0_wnn) - inverse_PH(chi_wnn)    
    return gamma_wnn

# ----------------------------------------------------------------------
def fixed_fermionic_window_python_wnk(chi_wnk, nwf):

    r""" Helper routine to reduce the number of fermionic Matsubara 
    frequencies :math:`\nu` in a two frequency and one momenta dependent
    generalized susceptibility :math:`\chi_{abcd}(\omega, \nu, \mathbf{k})`.

    Parameters
    ----------

    chi_wnk : two frequency and one momenta dependent generalized 
              susceptibility :math:`\chi_{abcd}(\omega, \nu, \mathbf{k})`.
    nwf : number of fermionic frequencies to keep.

    Returns
    -------

    chi_wnk_out : Susceptibility with reduced number of fermionic Matsubara
                  frequencies.
    """

    g2 = chi_wnk
    wmesh, nmesh, kmesh = g2.mesh.components
    
    beta = g2.mesh.components[0].beta
    nmesh_small = MeshImFreq(beta=beta, S='Fermion', n_max=nwf)

    chi_wnk_out = Gf(mesh=MeshProduct(wmesh, nmesh_small, kmesh), target_shape=g2.target_shape)

    n = g2.data.shape[1]
    s = n/2 - nwf
    e = n/2 + nwf
    
    chi_wnk_out.data[:] = g2.data[:, s:e, :]

    return chi_wnk_out

# ----------------------------------------------------------------------
def get_chi0_wnk(g_wk, nw=1, nwf=None):

    r""" Compute the generalized lattice bubble susceptibility 
    :math:`\chi^{(0)}_{abcd}(\omega, \nu, \mathbf{k})` from the single-particle
    Green's function :math:`G_{ab}(\omega, \mathbf{k})`.

    Parameters
    ----------

    g_wk : Single-particle Green's function :math:`G_{ab}(\omega, \mathbf{k})`.
    nw : Number of bosonic freqiencies in :math:`\chi`.
    nwf : Number of fermionic freqiencies in :math:`\chi`.    

    Returns
    -------

    chi0_wnk : Generalized lattice bubble susceptibility
               :math:`\chi^{(0)}_{abcd}(\omega, \nu, \mathbf{k})`
    """

    fmesh = g_wk.mesh.components[0]
    kmesh = g_wk.mesh.components[1]

    if nwf is None:
        nwf = len(fmesh) / 2

    mpi.barrier()
    mpi.report('g_wk ' + str(g_wk[Idx(2), Idx(0,1,2)][0,0]))
    n = np.sum(g_wk.data) / len(kmesh)
    mpi.report('n ' + str(n))
    mpi.barrier()

    mpi.report('--> g_wr from g_wk')
    g_wr = fourier_wk_to_wr(g_wk)

    mpi.barrier()
    mpi.report('g_wr ' + str(g_wr[Idx(2), Idx(0,1,2)][0,0]))
    n_r = np.sum(g_wr.data, axis=0)[0]
    mpi.report('n_r=0 ' + str(n_r[0,0]))
    mpi.barrier()
    
    mpi.report('--> chi0_wnr from g_wr')
    chi0_wnr = chi0r_from_gr_PH(nw=nw, nn=nwf, g_nr=g_wr)

    #mpi.report('--> chi0_wnr from g_wr (nompi)')
    #chi0_wnr_nompi = chi0r_from_gr_PH_nompi(nw=nw, nn=nwf, g_wr=g_wr)
    
    del g_wr

    #abs_diff = np.abs(chi0_wnr.data - chi0_wnr_nompi.data)
    #mpi.report('shape = ' + str(abs_diff.shape))
    #idx = np.argmax(abs_diff)
    #mpi.report('argmax = ' + str(idx))
    #diff = np.max(abs_diff)
    #mpi.report('diff = %6.6f' % diff)
    #del chi0_wnr
    #chi0_wnr = chi0_wnr_nompi

    #exit()
    
    mpi.barrier()
    mpi.report('chi0_wnr ' + str(chi0_wnr[Idx(0), Idx(0), Idx(0,0,0)][0,0,0,0]))
    chi0_r0 = np.sum(chi0_wnr[:, :, Idx(0,0,0)].data)
    mpi.report('chi0_r0 ' + str(chi0_r0))    
    mpi.barrier()
    
    mpi.report('--> chi0_wnk from chi0_wnr')
    chi0_wnk = chi0q_from_chi0r(chi0_wnr)

    del chi0_wnr

    mpi.barrier()
    mpi.report('chi0_wnk ' + str(chi0_wnk[Idx(0), Idx(0), Idx(0,0,0)][0,0,0,0]))
    chi0 = np.sum(chi0_wnk.data) / len(kmesh)
    mpi.report('chi0 = ' + str(chi0))
    mpi.barrier()

    #if mpi.is_master_node():
    if False:
        from triqs_tprf.ParameterCollection import ParameterCollection
        p = ParameterCollection()
        p.g_wk = g_wk
        p.g_wr = g_wr
        p.chi0_wnr = chi0_wnr
        p.chi0_wnk = chi0_wnk

        print '--> Writing debug info for BSE'
        with HDFArchive('data_debug_bse.h5', 'w') as arch:
            arch['p'] = p

    mpi.barrier()
    
    return chi0_wnk    
        
# ----------------------------------------------------------------------
def solve_lattice_bse(g_wk, gamma_wnn, tail_corr_nwf=None):

    r""" Compute the generalized lattice susceptibility 
    :math:`\chi_{abcd}(\omega, \mathbf{k})` using the Bethe-Salpeter 
    equation (BSE).

    Parameters
    ----------

    g_wk : Single-particle Green's function :math:`G_{ab}(\omega, \mathbf{k})`.
    gamma_wnn : Local particle-hole vertex function 
                :math:`\Gamma_{abcd}(\omega, \nu, \nu')`
    tail_corr_nwf : Number of fermionic freqiencies to use in the 
                    tail correction of the sum over fermionic frequencies.

    Returns
    -------

    chi0_wk : Generalized lattice susceptibility
              :math:`\chi_{abcd}(\omega, \mathbf{k})`
    """
    
    fmesh_g = g_wk.mesh.components[0]
    kmesh = g_wk.mesh.components[1]
    
    bmesh = gamma_wnn.mesh.components[0]
    fmesh = gamma_wnn.mesh.components[1]

    nk = len(kmesh)
    nw = (len(bmesh) + 1) / 2
    nwf = len(fmesh) / 2
    nwf_g = len(fmesh_g) / 2

    if mpi.is_master_node():
        print tprf_banner(), "\n"
        print 'Lattcie BSE with local vertex approximation.\n'
        print 'nk    =', nk
        print 'nw    =', nw
        print 'nwf   =', nwf
        print 'nwf_g =', nwf_g
        print 'tail_corr_nwf =', tail_corr_nwf
        print    

    if tail_corr_nwf is None:
        tail_corr_nwf = nwf

    mpi.report('--> chi0_wnk_tail_corr')
    chi0_wnk_tail_corr = get_chi0_wnk(g_wk, nw=nw, nwf=tail_corr_nwf)

    mpi.report('--> trace chi0_wnk_tail_corr (WARNING! NO TAIL FIT. FIXME!)')
    chi0_wk_tail_corr = chi0q_sum_nu_tail_corr_PH(chi0_wnk_tail_corr)
    #chi0_wk_tail_corr = chi0q_sum_nu(chi0_wnk_tail_corr)

    mpi.barrier()
    mpi.report('B1 ' + str(chi0_wk_tail_corr[Idx(0), Idx(0,0,0)][0,0,0,0]))    
    mpi.barrier()

    mpi.report('--> chi0_wnk_tail_corr to chi0_wnk')
    if tail_corr_nwf != nwf:
        mpi.report('--> fixed_fermionic_window_python_wnk')
        chi0_wnk = fixed_fermionic_window_python_wnk(chi0_wnk_tail_corr, nwf=nwf)
    else:
        chi0_wnk = chi0_wnk_tail_corr.copy()

    del chi0_wnk_tail_corr

    mpi.barrier()
    mpi.report('C ' + str(chi0_wnk[Idx(0), Idx(0), Idx(0,0,0)][0,0,0,0]))    
    mpi.barrier()
    
    mpi.report('--> trace chi0_wnk')
    chi0_wk = chi0q_sum_nu(chi0_wnk)

    mpi.barrier()
    mpi.report('D ' + str(chi0_wk[Idx(0), Idx(0,0,0)][0,0,0,0]))    
    mpi.barrier()

    dchi_wk = chi0_wk_tail_corr - chi0_wk

    chi0_kw = Gf(mesh=MeshProduct(kmesh, bmesh), target_shape=chi0_wk_tail_corr.target_shape)
    chi0_kw.data[:] = chi0_wk_tail_corr.data.swapaxes(0, 1)

    del chi0_wk
    del chi0_wk_tail_corr

    assert( chi0_wnk.mesh.components[0] == bmesh )
    assert( chi0_wnk.mesh.components[1] == fmesh )
    assert( chi0_wnk.mesh.components[2] == kmesh )

    # -- Lattice BSE calc with built in trace
    mpi.report('--> chi_kw from BSE')
    #mpi.report('DEBUG BSE INACTIVE'*72)
    chi_kw = chiq_sum_nu_from_chi0q_and_gamma_PH(chi0_wnk, gamma_wnn)
    #chi_kw = chi0_kw.copy()

    mpi.barrier()
    mpi.report('--> chi_kw from BSE (done)')
    
    del chi0_wnk

    mpi.report('--> chi_kw tail corrected (using chi0_wnk)')
    for k in kmesh:
        chi_kw[k, :] += dchi_wk[:, k] # -- account for high freq of chi_0 (better than nothing)

    del dchi_wk

    mpi.report('--> solve_lattice_bse, done.')

    return chi_kw, chi0_kw
 
# ----------------------------------------------------------------------
def solve_lattice_bse_depr(g_wk, gamma_wnn, tail_corr_nwf=-1):

    fmesh_huge, kmesh = g_wk.mesh.components
    bmesh = gamma_wnn.mesh.components[0]
    fmesh = gamma_wnn.mesh.components[1]

    nk = len(kmesh)
    nw = (len(bmesh) + 1) / 2
    nwf = len(fmesh) / 2
    nwf_sigma = len(fmesh_huge) / 2

    if mpi.is_master_node():    
        print tprf_banner(), "\n"

        print 'Lattcie BSE with local vertex approximation.\n'
        print 'nk  = ', nk
        print 'nw  = ', nw
        print 'nwf = ', nwf
        print 'nwf_sigma = ', nwf_sigma
        print 'nwf = ', tail_corr_nwf, ' (gf)'
        print    

    # -- Lattice BSE calc with built in trace using g_wk
    from triqs_tprf.lattice import chiq_sum_nu_from_g_wk_and_gamma_PH

    chi_kw = chiq_sum_nu_from_g_wk_and_gamma_PH(g_wk, gamma_wnn, tail_corr_nwf=tail_corr_nwf)
    
    return chi_kw
    
# ----------------------------------------------------------------------
def solve_lattice_bse_e_k_sigma_w(mu, e_k, sigma_w, gamma_wnn, tail_corr_nwf=-1):

    kmesh = e_k.mesh
    fmesh_huge = sigma_w.mesh
    bmesh = gamma_wnn.mesh.components[0]
    fmesh = gamma_wnn.mesh.components[1]

    nk = len(kmesh)
    nw = (len(bmesh) + 1) / 2
    nwf = len(fmesh) / 2
    nwf_sigma = len(fmesh_huge) / 2

    if mpi.is_master_node():    
        print tprf_banner(), "\n"

        print 'Lattcie BSE with local vertex approximation.\n'
        print 'nk  =', nk
        print 'nw  =', nw
        print 'nwf           =', nwf
        print 'nwf_sigma     =', nwf_sigma
        print 'nwf_chi0_tail =', tail_corr_nwf
        print    

    # -- Lattice BSE calc with built in trace using g_wk
    from triqs_tprf.lattice import chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH

    chi_kw = chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH(mu, e_k, sigma_w, gamma_wnn, tail_corr_nwf=tail_corr_nwf)
    
    return chi_kw
    
