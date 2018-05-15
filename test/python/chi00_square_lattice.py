
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.lattice import BrillouinZone, BravaisLattice
from pytriqs.lattice.tight_binding import TBLattice

from pytriqs.gf import MeshBrillouinZone
from pytriqs.gf import MeshCyclicLattice
from pytriqs.gf import MeshProduct
from pytriqs.gf import Gf, MeshImFreq, Idx

# ----------------------------------------------------------------------

from triqs_tprf.lattice import g0k_from_ek
from triqs_tprf.lattice import gr_from_gk
from triqs_tprf.lattice_utils import ek_tb_dispersion_on_bzmesh

# ----------------------------------------------------------------------

from triqs_tprf.lattice import chi00_wk_from_ek

# ----------------------------------------------------------------------

from triqs_tprf.lattice import grt_from_grw
from triqs_tprf.lattice import chi0_tr_from_grt_PH
from triqs_tprf.lattice import chi0_w0r_from_grt_PH
from triqs_tprf.lattice import chi_w0r_from_chi_tr
from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_wk_from_chi_wr

from triqs_tprf.lattice_utils import chi_w0r_from_chi_tr_np_trapz

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
    nnu = 400
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
    
    periodization_matrix = np.diag(np.array(list(n_k), dtype=np.int32))
    print 'periodization_matrix =\n', periodization_matrix

    print '--> tight binding model'
    tb_lattice = TBLattice(
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

    bz = BrillouinZone(tb_lattice.bl)
    bzmesh = MeshBrillouinZone(bz, periodization_matrix)
    ek = ek_tb_dispersion_on_bzmesh(tb_lattice, bzmesh, bz)
    
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw_g)
    lmesh = MeshCyclicLattice(tb_lattice.bl, periodization_matrix)

    print '--> g0k'
    g0_wk = g0k_from_ek(mu=mu, ek=ek, mesh=wmesh)

    print '--> g0r'
    g0_wr = gr_from_gk(g0_wk)

    # ------------------------------------------------------------------
    # -- anaytic chi00
    
    print '--> chi00wq analytic'
    chi00wq_analytic = chi00_wk_from_ek(ek, nw, beta, mu)

    # ------------------------------------------------------------------
    # -- imtime chi00

    print '--> grt from grw' 
    g0_tr = grt_from_grw(g0_wr)

    print '--> chi0_tr_from_grt_PH'
    chi00_tr = chi0_tr_from_grt_PH(g0_tr)

    print '--> chi0_tr_from_grt_PH'
    chi00_wr_opt = chi0_w0r_from_grt_PH(g0_tr)
    
    print '--> chi_wr_from_chi_tr'
    chi00_wr = chi_wr_from_chi_tr(chi00_tr, nw=1)

    print '--> chi_w0r_from_chi_tr'
    chi00_wr_ref1 = chi_w0r_from_chi_tr(chi00_tr)

    print '--> chi_wr_from_chi_tr using np.trapz'
    chi00_wr_ref2 = chi_w0r_from_chi_tr_np_trapz(chi00_tr)

    def comp(chi):
        return np.squeeze(chi[:, Idx(0, 0, 0)].data).reshape((4, 4))
    
    #print comp(chi00_wr).real
    #print comp(chi00_wr_ref1).real
    #print comp(chi00_wr_ref1).real
    
    #print comp(chi00_wr).real - comp(chi00_wr_ref1).real
    
    #print np.max(np.abs(chi00_wr.data - chi00_wr_ref1.data))

    np.testing.assert_array_almost_equal(chi00_wr.data, chi00_wr_opt.data, decimal=4)
    np.testing.assert_array_almost_equal(chi00_wr.data, chi00_wr_ref1.data, decimal=4)
    np.testing.assert_array_almost_equal(chi00_wr.data, chi00_wr_ref2.data, decimal=4)
    
    print '--> chi_wk_from_chi_wr'
    chi00wq_imtime = chi_wk_from_chi_wr(chi00_wr)
    
    # ------------------------------------------------------------------
    # -- imfreq chi00
    
    print '--> chi00r'
    chi00r = chi0r_from_gr_PH(nw=nw, nnu=nnu, gr=g0_wr)

    print '--> chi00q'
    chi00q = chi0q_from_chi0r(chi00r)

    print '--> chi00wq'
    chi00wq_imfreq = chi0q_sum_nu(chi00q)

    # ------------------------------------------------------------------
    # -- Compare results

    def cf_chi_w0(chi1, chi2, decimal=9):
        chi1, chi2 = chi1[Idx(0), :].data, chi2[Idx(0), :].data
        diff = np.linalg.norm(chi1 - chi2)
        print '|dchi| =', diff
        np.testing.assert_array_almost_equal(chi1, chi2, decimal=decimal)
    
    print '--> Cf analytic with imtime'
    cf_chi_w0(chi00wq_analytic, chi00wq_imtime, decimal=4)

    print '--> Cf analytic with imfreq'
    cf_chi_w0(chi00wq_analytic, chi00wq_imfreq, decimal=2)
    
# ----------------------------------------------------------------------
if __name__ == '__main__':

    test_square_lattice_chi00()
