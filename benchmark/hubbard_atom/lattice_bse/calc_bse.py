# ----------------------------------------------------------------------

import glob
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import *
from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt     
from pytriqs.lattice import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection

from triqs_tprf.analytic_hubbard_atom import analytic_hubbard_atom

# ----------------------------------------------------------------------

from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued
from triqs_tprf.freq_conv import from_3nu_PH

from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.chi_from_gg2 import chi_from_gg2_PH

from triqs_tprf.linalg import inverse_PH

# ----------------------------------------------------------------------

from triqs_tprf.lattice import g0k_from_ek
from triqs_tprf.lattice import gk_from_ek_sigma
from triqs_tprf.lattice import gr_from_gk

from triqs_tprf.lattice import chi0r_from_gr_PH
from triqs_tprf.lattice import chi0q_from_chi0r
from triqs_tprf.lattice import chi0q_sum_nu
from triqs_tprf.lattice import chi0q_sum_nu_q

from triqs_tprf.lattice import chiq_from_chi0q_and_gamma_PH
from triqs_tprf.lattice import chiq_sum_nu, chiq_sum_nu_q

# ----------------------------------------------------------------------
def setup_dispersion_and_bz(parm):

    """Ugly hack for evaluating the dispersion of a tight binding lattice
    on a MeshBrillouinZone Green's function mesh."""

    nk = 1
    print 'nk =', nk
    
    # -- Get the bravais lattice from parm and constr the Brillouin zone

    bz = BrillouinZone(BravaisLattice([(1,0,0), (0,1,0)]))
    
    bzmesh = MeshBrillouinZone(bz, n_k=nk)
    q_list = np.array([q.value[:2] for q in bzmesh])
    q_list /= 2.*np.pi

    #for qidx, q in enumerate(bzmesh):
    #    print qidx, q.value/(2.*np.pi)

    # ------------------------------------------------------------------
    # -- Setup dispersion

    # using zero dispersion, to get atomic limit
    
    ek = Gf(mesh=bzmesh, target_shape=[1, 1])

    for kidx, k in enumerate(bzmesh):
        k_val = k.value[:, None] / (2. * np.pi)
        kidx_0 = 0
        ek.data[kidx, :, :] = 0.0

    nk_tot = len(bzmesh)
    print 'nk_tot =', nk_tot
    
    return ek, bz, nk_tot
        
# ----------------------------------------------------------------------
def solve_lattice_bse(parm, momsus=False):

    print '--> solve_lattice_bse'

    print 'nw =', parm.nw
    print 'nwf =', parm.nwf
    
    # ------------------------------------------------------------------
    # -- Setup lattice

    ek, bz, parm.nk = setup_dispersion_and_bz(parm)
          
    # ------------------------------------------------------------------
    # -- Lattice single-particle Green's function

    mesh = MeshImFreq(beta=parm.beta, S='Fermion', n_max=parm.nwf_gf)

    parm.Sigma_iw = parm.G_iw.copy()
    G0_iw = parm.G_iw.copy()

    G0_iw << inverse(iOmega_n + 0.5*parm.U)
    parm.Sigma_iw << inverse(G0_iw) - inverse(parm.G_iw)
    
    parm.mu = 0.5*parm.U
    gk = gk_from_ek_sigma(mu=parm.mu, ek=ek, sigma=parm.Sigma_iw)
    gr = gr_from_gk(gk)    

    # ------------------------------------------------------------------
    # -- check that G_imp == G_loc

    # -- Get the 0,0,0 real space component of the Green's function

    if False:
        # Complicated way
        rmesh = gr.mesh.components[1]
        r_mesh_points = np.array([r for r in rmesh])
        ridx = rmesh.index_to_linear((0,0,0))
        r_mesh_point = r_mesh_points[ridx]
        G_iw = gr[:, r_mesh_point]
    else:
        # (Less) Complicated way, thanks Nils Wentzell! :)
        rmesh = gr.mesh.components[1]
        ridx = rmesh.index_to_linear((0,0,0))
        G_iw = gr[:, MeshPoint(ridx)]
    
    #oplot(parm.G_iw)
    #oplot(G_iw)    
    #plt.show()

    if False:
        diff = G_iw.data - parm.G_iw.data
        print '--> Checking G_imp == G_loc,',
        print 'max|diff| =', np.max(np.abs(diff))

    np.testing.assert_array_almost_equal(G_iw.data, parm.G_iw.data)

    # ------------------------------------------------------------------
    # -- Non-interacting generalized lattice susceptibility

    chi0r = chi0r_from_gr_PH(nw=parm.nw, nnu=parm.nwf, gr=gr) 
    chi0q = chi0q_from_chi0r(chi0r, bz)

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

    beta = 1.234
    nwf = 248
    
    p = ParameterCollection(
        beta = beta,
        U = 5.0,
        nw = 1,
        nwf = nwf,
        nwf_gf = 2*nwf,
        )

    p = analytic_hubbard_atom(**p.dict())
    p = solve_lattice_bse(p)

    # -- Tracing the analytic generalized susceptibility
    p.chi_w_analytic = np.sum(p.chi_m.data) / p.beta**2    

    # -- Tracing the numerical generalized susceptibility
    p.chi_w_ref = np.sum(p.chiq.data) / p.beta**2 / p.nk

    print 'chi_w_analytic =', p.chi_w_analytic
    print 'chi_w.data     =', p.chi_w.data.flatten()
    print 'chi_w_ref      =', p.chi_w_ref
    print 'chi_m_static   =', p.chi_m_static
    
    np.testing.assert_array_almost_equal(
        np.squeeze(p.chi_w.data), p.chi_w_ref)
    
