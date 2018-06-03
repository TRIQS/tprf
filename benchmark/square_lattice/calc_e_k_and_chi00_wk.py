
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lindhard_chi00_wk

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    n_k = (64, 64, 1)
    nw = 100
    
    beta = 5.0
    mu = 0.0
    t = 1.0
    
    print '--> tight binding model'
    T = - t * np.eye(1)
    t_r = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        hopping = {
            # nearest neighbour hopping -t
            ( 0,+1): T,
            ( 0,-1): T,
            (+1, 0): T,
            (-1, 0): T,
            },
        orbital_positions = [(0,0,0)],
        orbital_names = ['0'],
        )

    print '--> dispersion e_k'
    e_k = t_r.on_mesh_brillouin_zone(n_k)

    print '--> lattice g0_wk'
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

    # -- Call TPRF chi0_wk bubble calc
    chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)

    # -- Analytic reference (Lindhard formula in momentum space)
    print '--> lindhard_chi00_wk'
    n_k = (32, 32, 1)
    e_k = t_r.on_mesh_brillouin_zone(n_k)
    chi00_wk_analytic = lindhard_chi00_wk(e_k=e_k, nw=1, beta=beta, mu=mu)

    filename = 'data_e_k_and_chi00_wk.h5'
    
    with HDFArchive(filename, 'w') as arch:
        arch['e_k'] = e_k
        arch['chi00_wk'] = chi00_wk
        arch['chi00_wk_analytic'] = chi00_wk_analytic
