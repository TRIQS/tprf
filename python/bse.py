
import numpy as np

# ----------------------------------------------------------------------

import pytriqs.utility.mpi as mpi
from pytriqs.gf import MeshImFreq, MeshProduct, Gf, Idx

# ----------------------------------------------------------------------

from triqs_tprf.logo import tprf_banner

from triqs_tprf.linalg import inverse_PH

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi0r_from_gr_PH
from triqs_tprf.lattice import chi0q_from_chi0r
from triqs_tprf.lattice import chi0q_sum_nu, chi0q_sum_nu_tail_corr_PH
from triqs_tprf.lattice import chiq_sum_nu_from_chi0q_and_gamma_PH

# ----------------------------------------------------------------------
def solve_local_bse(chi0_wnn, chi_wnn):

    gamma_wnn = inverse_PH(chi0_wnn) - inverse_PH(chi_wnn)
    
    return gamma_wnn

# ----------------------------------------------------------------------
def fixed_fermionic_window_python_wnk(g2, nwf):

    #nw = (g2.data.shape[0] + 1) / 2

    wmesh, nmesh, kmesh = g2.mesh.components
    
    beta = g2.mesh.components[0].beta
    nmesh_small = MeshImFreq(beta=beta, S='Fermion', n_max=nwf)

    g2_out = Gf(mesh=MeshProduct(wmesh, nmesh_small, kmesh), target_shape=g2.target_shape)

    n = g2.data.shape[1]
    s = n/2 - nwf
    e = n/2 + nwf
    
    g2_out.data[:] = g2.data[:, s:e, :]

    return g2_out

# ----------------------------------------------------------------------
def get_chi0_wnk(g_wk, nw=1, nwf=None):

    fmesh = g_wk.mesh.components[0]

    if nwf is None:
        nwf = len(fmesh) / 2

    mpi.report('--> g_wr from g_wk')
    g_wr = fourier_wk_to_wr(g_wk)

    mpi.report('--> chi0_wnr from g_wr')
    chi0_wnr = chi0r_from_gr_PH(nw=nw, nnu=nwf, gr=g_wr)

    del g_wr

    mpi.report('--> chi0_wnk from chi0_wnr')
    chi0_wnk = chi0q_from_chi0r(chi0_wnr)

    del chi0_wnr

    return chi0_wnk    
        
# ----------------------------------------------------------------------
def solve_lattice_bse(g_wk, gamma_wnn, tail_corr_nwf=None):

    kmesh = g_wk.mesh.components[1]
    bmesh = gamma_wnn.mesh.components[0]
    fmesh = gamma_wnn.mesh.components[1]

    nk = len(kmesh)
    nw = (len(bmesh) + 1) / 2
    nwf = len(fmesh) / 2

    if mpi.is_master_node():
        print tprf_banner(), "\n"
        print 'Lattcie BSE with local vertex approximation.\n'
        print 'nk  = ', nk
        print 'nw  = ', nw
        print 'nwf = ', nwf
        print    

    if tail_corr_nwf is None:
        tail_corr_nwf = nwf

    mpi.report('--> chi0_wnk_tail_corr')
    chi0_wnk_tail_corr = get_chi0_wnk(g_wk, nw=nw, nwf=tail_corr_nwf)

    mpi.report('--> trace chi0_wnk_tail_corr')
    chi0_wk_tail_corr = chi0q_sum_nu_tail_corr_PH(chi0_wnk_tail_corr)

    mpi.report('--> chi0_wnk_tail_corr to chi0_wnk')
    if tail_corr_nwf != nwf:
        mpi.report('--> fixed_fermionic_window_python_wnk')
        chi0_wnk = fixed_fermionic_window_python_wnk(chi0_wnk_tail_corr, nwf=nwf)
    else:
        chi0_wnk = chi0_wnk_tail_corr.copy()

    del chi0_wnk_tail_corr

    mpi.report('--> trace chi0_wnk')
    chi0_wk = chi0q_sum_nu(chi0_wnk)

    dchi_wk = chi0_wk_tail_corr - chi0_wk

    del chi0_wk
    del chi0_wk_tail_corr

    assert( chi0_wnk.mesh.components[0] == bmesh )
    assert( chi0_wnk.mesh.components[1] == fmesh )
    assert( chi0_wnk.mesh.components[2] == kmesh )

    # -- Lattice BSE calc with built in trace
    mpi.report('--> chi_kw from BSE')
    chi_kw = chiq_sum_nu_from_chi0q_and_gamma_PH(chi0_wnk, gamma_wnn)

    del chi0_wnk

    mpi.report('--> chi_kw tail corrected (using chi0_wnk)')
    for k in chi_kw.mesh.components[0]:
        chi_kw[k, :] += dchi_wk[:, k] # -- account for high freq of chi_0 (better than nothing)

    del dchi_wk

    mpi.report('--> solve_lattice_bse, done.')

    return chi_kw
 
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
    
