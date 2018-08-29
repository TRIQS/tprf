
import numpy as np

# ----------------------------------------------------------------------

import pytriqs.utility.mpi as mpi

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
def get_chi0_wnk(g_wk, nw=1, nwf=None):

    fmesh = g_wk.mesh.components[0]

    if nwf is None:
        nwf = len(fmesh) / 2

    print '--> g_wr from g_wk'
    g_wr = fourier_wk_to_wr(g_wk)

    print '--> chi0_wnr from g_wr'
    chi0_wnr = chi0r_from_gr_PH(nw=nw, nnu=nwf, gr=g_wr)

    del g_wr

    print '--> chi0_wnk from chi0_wnr'
    chi0_wnk = chi0q_from_chi0r(chi0_wnr)

    return chi0_wnk    
        
# ----------------------------------------------------------------------
def solve_lattice_bse_depr(g_wk, gamma_wnn, tail_corr_nwf=None, return_chi0_wk=False):

    kmesh = g_wk.mesh.components[1]
    bmesh = gamma_wnn.mesh.components[0]
    fmesh = gamma_wnn.mesh.components[1]

    nk = len(kmesh)
    nw = (len(bmesh) + 1) / 2
    nwf = len(fmesh) / 2

    print tprf_banner(), "\n"

    print 'Lattcie BSE with local vertex approximation.\n'
    print 'nk  = ', nk
    print 'nw  = ', nw
    print 'nwf = ', nwf
    print    

    chi0_wnk = get_chi0_wnk(g_wk, nw=nw, nwf=nwf)
    print '--> trace chi0_wnk'
    chi0_wk = chi0q_sum_nu(chi0_wnk)
    
    assert( chi0_wnk.mesh.components[0] == bmesh )
    assert( chi0_wnk.mesh.components[1] == fmesh )
    assert( chi0_wnk.mesh.components[2] == kmesh )

    # -- Trace chi0_kw to get correction to chi_kw

    print '--> trace chi0_wnk (tail corr)'
    if tail_corr_nwf is None:
        chi0_wk_tail_corr = chi0q_sum_nu_tail_corr_PH(chi0_wnk)
    else:
        chi0_wnk_tail_corr = get_chi0_wnk(g_wk, nw=nw, nwf=tail_corr_nwf)
        chi0_wk_tail_corr = chi0q_sum_nu_tail_corr_PH(chi0_wnk_tail_corr)

    dchi_wk = chi0_wk_tail_corr - chi0_wk

    del chi0_wk
    #del chi0_wk_tail_corr

    if tail_corr_nwf is not None:
        del chi0_wnk_tail_corr
    
    # -- Lattice BSE calc with built in trace
    print '--> chi_kw from BSE'
    chi_kw = chiq_sum_nu_from_chi0q_and_gamma_PH(chi0_wnk, gamma_wnn)

    del chi0_wnk

    print '--> chi_kw tail corrected (using chi0_wnk)'
    for k in chi_kw.mesh.components[0]:
        chi_kw[k, :] += dchi_wk[:, k] # -- account for high freq of chi_0 (better than nothing)

    del dchi_wk

    if return_chi0_wk:
        return chi_kw, chi0_wk_tail_corr
    else:
        return chi_kw

# ----------------------------------------------------------------------
def solve_lattice_bse(g_wk, gamma_wnn, tail_corr_nwf=-1):

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
    
