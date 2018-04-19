
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive

from pytriqs.lattice import BrillouinZone, BravaisLattice
from pytriqs.lattice.tight_binding import TBLattice

from pytriqs.gf import Gf
from pytriqs.gf import MeshImFreq
from pytriqs.gf import MeshProduct
from pytriqs.gf import MeshBrillouinZone
from pytriqs.gf import MeshCyclicLattice

# ----------------------------------------------------------------------

from triqs_tprf.lattice import chi00_wk_from_ek
from triqs_tprf.lattice_utils import ek_tb_dispersion_on_bzmesh

from triqs_tprf.lattice import g0k_from_ek
from triqs_tprf.lattice import gr_from_gk
from triqs_tprf.lattice import grt_from_grw
from triqs_tprf.lattice import chi0_tr_from_grt_PH
from triqs_tprf.lattice import chi0_w0r_from_grt_PH
from triqs_tprf.lattice import chi_wr_from_chi_tr
from triqs_tprf.lattice import chi_wk_from_chi_wr

from triqs_tprf.lattice_utils import chi_w0r_from_chi_tr_np_trapz
    
# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    #n_k = (128, 128, 1) # need some parallelisation for this... FIXME
    #n_k = (96, 96, 1)
    n_k = (64, 64, 1)
    #n_k = (32, 32, 1)
    #n_k = (16, 16, 1)
    #n_k = (6, 6, 1)
    #n_k = (2, 2, 1)
    nw = 1
    
    # ------------------------------------------------------------------
    # -- tight binding parameters

    beta = 20.0
    mu = 0.0
    t = 1.0
    
    h_loc = np.zeros((2, 2))        
    T = - t * np.eye(2)

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

    n_orb = tb_lattice.NOrbitalsInUnitCell
    bz = BrillouinZone(tb_lattice.bl)
    bzmesh = MeshBrillouinZone(bz, periodization_matrix)
    ek = ek_tb_dispersion_on_bzmesh(tb_lattice, bzmesh, bz)
    
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=100)

    print '--> g0k'
    g0k = g0k_from_ek(mu=mu, ek=ek, mesh=wmesh)

    print '--> g0r'
    g0r = gr_from_gk(g0k)

    print '--> grt from grw' 
    g0rt = grt_from_grw(g0r)

    print '--> chi0_w0r_from_grt_PH (bubble in tau & r)'
    chi00_wr = chi0_w0r_from_grt_PH(g0rt)

    #print '--> chi0_tr_from_grt_PH'
    #chi00_tr = chi0_tr_from_grt_PH(g0rt)

    #print '--> chi_wr_from_chi_tr'
    #chi00_wr = chi_wr_from_chi_tr(chi00_tr, nw=1) # this could be replaced by integral over tau only.

    #print '--> chi_wr_from_chi_tr using np.trapz'
    #chi00_wr = chi_w0r_from_chi_tr_np_trapz(chi00_tr)
    
    print '--> chi_wk_from_chi_wr'
    chi00wq_imtime = chi_wk_from_chi_wr(chi00_wr)

    # -- Analytic form
    
    print '--> chi00wq'
    chi00wq = chi00_wk_from_ek(ek, nw, beta, mu)

    filename = 'data_ek_and_chi00wq.h5'
    
    with HDFArchive(filename, 'w') as arch:
        arch['ek'] = ek
        arch['chi00wq'] = chi00wq
        arch['chi00wq_imtime'] = chi00wq_imtime
