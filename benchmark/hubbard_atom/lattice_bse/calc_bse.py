# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import Gf
from pytriqs.gf import inverse
from pytriqs.gf import iOmega_n
from pytriqs.gf import MeshImFreq
from pytriqs.gf import MeshBrillouinZone

from pytriqs.lattice import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

#from triqs_tprf.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------

from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.lattice import fourier_wk_to_wr

from triqs_tprf.lattice import chi0r_from_gr_PH
from triqs_tprf.lattice import chi0q_from_chi0r
from triqs_tprf.lattice import chi0q_sum_nu
from triqs_tprf.lattice import chi0q_sum_nu_q

from triqs_tprf.lattice import chiq_from_chi0q_and_gamma_PH
from triqs_tprf.lattice import chiq_sum_nu, chiq_sum_nu_q

# ----------------------------------------------------------------------
        
from triqs_tprf.analytic_hubbard_atom import analytic_hubbard_atom

# ----------------------------------------------------------------------
def solve_lattice_bse(parm, momsus=False):

    print '--> solve_lattice_bse'

    print 'nw =', parm.nw
    print 'nwf =', parm.nwf
    
    # ------------------------------------------------------------------
    # -- Setup lattice

    bl = BravaisLattice([(1,0,0), (0,1,0)])

    bz = BrillouinZone(bl)    
    bzmesh = MeshBrillouinZone(bz, n_k=1) # only one k-point
    
    e_k = Gf(mesh=bzmesh, target_shape=[1, 1])
    e_k *= 0.
          
    # ------------------------------------------------------------------
    # -- Lattice single-particle Green's function

    mesh = MeshImFreq(beta=parm.beta, S='Fermion', n_max=parm.nwf_gf)

    parm.Sigma_iw = parm.G_iw.copy()
    G0_iw = parm.G_iw.copy()

    G0_iw << inverse(iOmega_n + 0.5*parm.U)
    parm.Sigma_iw << inverse(G0_iw) - inverse(parm.G_iw)
    
    parm.mu = 0.5*parm.U
    g_wk = lattice_dyson_g_wk(mu=parm.mu, e_k=e_k, sigma_w=parm.Sigma_iw)
    g_wr = fourier_wk_to_wr(g_wk)

    # ------------------------------------------------------------------
    # -- Non-interacting generalized lattice susceptibility

    chi0r = chi0r_from_gr_PH(nw=parm.nw, nnu=parm.nwf, gr=gr) 
    chi0q = chi0q_from_chi0r(chi0r)

    # ------------------------------------------------------------------
    # -- Solve lattice BSE

    parm.chiq = chiq_from_chi0q_and_gamma_PH(chi0q, parm.gamma_m)

    # ------------------------------------------------------------------
    # -- Store results and static results
    
    num = np.squeeze(parm.chiq.data.real)
    ref = np.squeeze(parm.chi_m.data.real)
    
    diff = np.max(np.abs(num - ref))
    print 'diff =', diff
    
    parm.chi_w = chiq_sum_nu_q(parm.chiq) # static suscept
    
    return parm

# ----------------------------------------------------------------------
if __name__ == '__main__':

    p = analytic_hubbard_atom(
        beta = 1.234,
        U = 5.0,
        nw = 1,
        nwf = 248,
        nwf_gf = 2 * 248,
        )
    
    p = solve_lattice_bse(p)

    # -- Tracing the analytic generalized susceptibility
    p.chi_w_analytic = np.sum(p.chi_m.data) / p.beta**2    

    # -- Tracing the numerical generalized susceptibility
    p.chi_w_ref = np.sum(p.chiq.data) / p.beta**2

    print 'chi_w_analytic =', p.chi_w_analytic
    print 'chi_w.data     =', p.chi_w.data.flatten()
    print 'chi_w_ref      =', p.chi_w_ref
    print 'chi_m_static   =', p.chi_m_static
    
    np.testing.assert_array_almost_equal(
        np.squeeze(p.chi_w.data), p.chi_w_ref)
    
