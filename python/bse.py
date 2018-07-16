
import numpy as np

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
def solve_lattice_bse(g_wk, gamma_wnn, nw=1, nwf=None):

    print tprf_banner(), "\n"

    print 'Lattcie BSE with local vertex approximation.\n'
    print 'nk  = ', g_wk.data.shape[1]
    print 'nw  = ', nw
    print 'nwf = ', nwf
    print    
    
    print '--> g_wr from g_wk'
    g_wr = fourier_wk_to_wr(g_wk)

    if nwf is None:
        raise ValueError, "set the number of fermionic frequencies, nwf"
        
    print '--> chi0_wnr from g_wr'
    chi0_wnr = chi0r_from_gr_PH(nw=nw, nnu=nwf, gr=g_wr)

    del g_wr

    print '--> chi0_wnk from chi0_wnr'
    chi0_wnk = chi0q_from_chi0r(chi0_wnr)

    del chi0_wnr

    # -- Trace chi0_kw to get correction to chi_kw

    print '--> trace chi0_wnk'
    chi0_wk = chi0q_sum_nu(chi0_wnk)
    print '--> trace chi0_wnk (tail corr)'
    chi0_wk_tail_corr = chi0q_sum_nu_tail_corr_PH(chi0_wnk)

    dchi_wk = chi0_wk_tail_corr - chi0_wk

    del chi0_wk
    del chi0_wk_tail_corr
    
    # -- Lattice BSE calc with built in trace
    print '--> chi_kw from BSE'
    chi_kw = chiq_sum_nu_from_chi0q_and_gamma_PH(chi0_wnk, gamma_wnn)

    del chi0_wnk

    print '--> chi_kw tail corrected (using chi0_wnk)'
    for k in chi_kw.mesh.components[0]:
        chi_kw[k, :] += dchi_wk[:, k] # -- account for high freq of chi_0 (better than nothing)

    del dchi_wk
    
    return chi_kw
