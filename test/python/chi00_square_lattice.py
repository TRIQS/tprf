
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq, Idx

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g0_wk

# ----------------------------------------------------------------------

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr

from triqs_tprf.lattice import chi0_tr_from_grt_PH
from triqs_tprf.lattice import chi0_w0r_from_grt_PH
from triqs_tprf.lattice import chi_w0r_from_chi_tr
from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_wk_from_chi_wr
from triqs_tprf.lattice import chi_wr_from_chi_wk

from triqs_tprf.lattice import lindhard_chi00_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------

from triqs_tprf.lattice import chi0r_from_gr_PH
from triqs_tprf.lattice import chi0q_from_chi0r
from triqs_tprf.lattice import chi0q_sum_nu
from triqs_tprf.lattice import chi0q_sum_nu_tail_corr_PH

# ----------------------------------------------------------------------
def test_square_lattice_chi00():
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    n_k = (2, 2, 1)
    nw_g = 500
    nn = 400
    nw = 1
    
    # ------------------------------------------------------------------
    # -- tight binding parameters

    beta = 20.0
    mu = 0.0
    t = 1.0
    
    h_loc = np.array([
        [-0.3, -0.5],
        [-0.5, .4],
        ])
        
    T = - t * np.array([
        [1., 0.23],
        [0.23, 0.5],
        ])

    # ------------------------------------------------------------------
    # -- tight binding
    
    print('--> tight binding model')
    t_r = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        hopping = {
            # nearest neighbour hopping -t
            ( 0, 0): h_loc,
            ( 0,+1): T,
            ( 0,-1): T,
            (+1, 0): T,
            (-1, 0): T,
            },
        orbital_positions = [(0,0,0)]*2,
        orbital_names = ['up_0', 'do_0'],
        )

    kmesh = t_r.get_kmesh(n_k)
    e_k = t_r.fourier(kmesh)

    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw_g)

    print('--> g0_wk')
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

    print('--> g0_wr')
    g0_wr = fourier_wk_to_wr(g0_wk)

    print('--> g0_tr') 
    g0_tr = fourier_wr_to_tr(g0_wr)
    
    # ------------------------------------------------------------------
    # -- anaytic chi00
    
    print('--> chi00_wk analytic')
    chi00_wk_analytic = lindhard_chi00_wk(e_k=e_k, nw=nw, beta=beta, mu=mu)

    print('--> chi00_wr analytic')
    chi00_wr_analytic = chi_wr_from_chi_wk(chi00_wk_analytic)
    
    # ------------------------------------------------------------------
    # -- imtime chi00

    print('--> chi0_tr_from_grt_PH')
    chi00_tr = chi0_tr_from_grt_PH(g0_tr)
    
    print('--> chi_wr_from_chi_tr')
    chi00_wr = chi_wr_from_chi_tr(chi00_tr, nw=1)

    print('--> chi_w0r_from_chi_tr')
    chi00_wr_ref = chi_w0r_from_chi_tr(chi00_tr)

    print('--> chi0_w0r_from_grt_PH')
    chi00_wr_opt = chi0_w0r_from_grt_PH(g0_tr)

    print('dchi00_wr     =', np.max(np.abs(chi00_wr_analytic.data - chi00_wr.data)))
    print('dchi00_wr_ref =', np.max(np.abs(chi00_wr_analytic.data - chi00_wr_ref.data)))
    print('dchi00_wr_opt =', np.max(np.abs(chi00_wr_analytic.data - chi00_wr_opt.data)))

    np.testing.assert_array_almost_equal(
        chi00_wr_analytic.data, chi00_wr.data, decimal=8)

    np.testing.assert_array_almost_equal(
        chi00_wr_analytic.data, chi00_wr_ref.data, decimal=4)

    np.testing.assert_array_almost_equal(
        chi00_wr_analytic.data, chi00_wr_opt.data, decimal=4)

    print('--> chi_wk_from_chi_wr')
    chi00_wk_imtime = chi_wk_from_chi_wr(chi00_wr)
    
    # ------------------------------------------------------------------
    # -- imtime chi00 helper function

    chi00_wk_imtime_2 = imtime_bubble_chi0_wk(g0_wk, nw=1)    

    # ------------------------------------------------------------------
    # -- imfreq chi00
    
    print('--> chi00_wnr')
    chi00_wnr = chi0r_from_gr_PH(nw=1, nn=nn, g_nr=g0_wr)

    print('--> chi00_wnk')
    chi00_wnk = chi0q_from_chi0r(chi00_wnr)

    # -- Test per k and w calculator for chi0_wnk
    print('--> chi00_wnk_ref')
    from triqs_tprf.lattice import chi0q_from_g_wk_PH
    chi00_wnk_ref = chi0q_from_g_wk_PH(nw=1, nn=nn, g_wk=g0_wk)

    diff = np.max(np.abs(chi00_wnk.data - chi00_wnk_ref.data))
    print('chi00_wnk diff =', diff)
    np.testing.assert_array_almost_equal(chi00_wnk.data, chi00_wnk_ref.data)
    
    print('--> chi00_wk_imfreq')
    chi00_wk_imfreq = chi0q_sum_nu(chi00_wnk)

    print('--> chi00_wk_imfreq_tail_corr')
    chi00_wk_imfreq_tail_corr = chi0q_sum_nu_tail_corr_PH(chi00_wnk)

    # ------------------------------------------------------------------
    # -- Compare results

    def cf_chi_w0(chi1, chi2, decimal=9):
        chi1, chi2 = chi1[Idx(0), :].data, chi2[Idx(0), :].data
        diff = np.linalg.norm(chi1 - chi2)
        print('|dchi| =', diff)
        np.testing.assert_array_almost_equal(chi1, chi2, decimal=decimal)
    
    print('--> Cf analytic with imtime')
    cf_chi_w0(chi00_wk_analytic, chi00_wk_imtime, decimal=7)

    print('--> Cf analytic with imtime 2')
    cf_chi_w0(chi00_wk_analytic, chi00_wk_imtime_2, decimal=4)
    
    print('--> Cf analytic with imfreq')
    cf_chi_w0(chi00_wk_analytic, chi00_wk_imfreq, decimal=2)

    print('--> Cf analytic with imfreq (tail corr)')
    cf_chi_w0(chi00_wk_analytic, chi00_wk_imfreq_tail_corr, decimal=5)
    
# ----------------------------------------------------------------------
if __name__ == '__main__':

    test_square_lattice_chi00()
