
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive

from pytriqs.lattice import BrillouinZone, BravaisLattice
from pytriqs.lattice.tight_binding import TBLattice

from pytriqs.gf import MeshBrillouinZone
from pytriqs.gf import MeshCyclicLattice

# ----------------------------------------------------------------------

from triqs_tprf.lattice import chi00_wk_from_ek
from triqs_tprf.lattice_utils import ek_tb_dispersion_on_bzmesh

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    #n_k = (128, 128, 1) # need some parallelisation for this... FIXME
    n_k = (64, 64, 1)
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

    bz = BrillouinZone(tb_lattice.bl)
    bzmesh = MeshBrillouinZone(bz, periodization_matrix)
    ek = ek_tb_dispersion_on_bzmesh(tb_lattice, bzmesh, bz)
    
    print '--> chi00wq'
    chi00wq = chi00_wk_from_ek(ek, nw, beta, mu)

    filename = 'data_ek_and_chi00wq.h5'
    
    with HDFArchive(filename, 'w') as arch:
        arch['ek'] = ek
        arch['chi00wq'] = chi00wq
