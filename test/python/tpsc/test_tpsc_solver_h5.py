
import numpy as np

from h5 import HDFArchive

from triqs.gfs import Gf
from triqs.mesh import MeshDLRImFreq, MeshImFreq
from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.tpsc_solver import tpsc_solver

def test_tpsc_solver_h5():

    beta = 15.0
    U = 2.0
    t = 1.0
    n = 1.0/2
    
    wmesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=12.0, eps=1e-14)

    tb = TBLattice(
        units=[(1,0,0), (0,1,0)],
        hoppings={(+1,+0) : [[-t]],
                  (-1,+0) : [[-t]],
                  (+0,+1) : [[-t]],
                  (+0,-1) : [[-t]]})

    e_k = tb.fourier(tb.get_kmesh(n_k=8))

    S = tpsc_solver(n, U, wmesh, e_k)
    S.solve()

    filename = 'data_h5_io.h5'
    with HDFArchive(filename, 'w') as A: A['S'] = S
    with HDFArchive(filename, 'r') as A: S_ref = A['S']
          
    assert( S == S_ref )


if __name__ == '__main__':

    test_tpsc_solver_h5()
