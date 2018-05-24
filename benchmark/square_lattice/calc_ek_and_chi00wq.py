
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
    
    n_k = (32, 32, 1)
    nw = 100
    
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
    e_k = ek_tb_dispersion_on_bzmesh(tb_lattice, bzmesh, bz)
    
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)

    g0_wk = g0k_from_ek(mu=mu, ek=e_k, mesh=wmesh)
    g0_wr = gr_from_gk(g0_wk)
    g0_tr = grt_from_grw(g0_wr)

    chi00_wr = chi0_w0r_from_grt_PH(g0_tr)
    chi00_wk = chi_wk_from_chi_wr(chi00_wr)

    # -- Analytic form
    chi00_wk_analytic = chi00_wk_from_ek(ek_in=e_k, nw=1, beta=beta, mu=mu)

    filename = 'data_e_k_and_chi00_wk.h5'
    
    with HDFArchive(filename, 'w') as arch:
        arch['e_k'] = e_k
        arch['chi00_wk'] = chi00_wk
        arch['chi00_wk_analytic'] = chi00_wk_analytic
