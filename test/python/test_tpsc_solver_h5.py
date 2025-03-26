
import numpy as np

from h5 import HDFArchive

from triqs.gf import MeshDLRImFreq, MeshImFreq, Gf
from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.tpsc_solver import tpsc_solver

def test_tpsc_solver_h5():

    beta = 15.0
    wmesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=12.0, eps=1e-14)

    U = 2.0
    n = 1.0/2
    units = [(1,0,0), (0,1,0)]
    hoppings = {(+1,+0) : [[-1.0]],
                (-1,+0) : [[-1.0]],
                (+0,+1) : [[-1.0]],
                (+0,-1) : [[-1.0]]}

    Lat = TBLattice(units=units, hoppings=hoppings)
    kmesh = Lat.get_kmesh(n_k=16)
    e_k = Lat.fourier(kmesh)

    # initialize solver
    S = tpsc_solver(n, U, wmesh, e_k, verbose=True)

    # run solver
    S.solve(calc_sigma=False, calc_g=False, check_self_consistency=False)

    Usp1 = S.Usp
    Uch1 = S.Uch
    docc1 = S.docc

    filename = 'data_h5_io.h5'

    with HDFArchive(filename, 'w') as A: A['S'] = S
    with HDFArchive(filename, 'r') as A: S_ref = A['S']
          
    assert( S == S_ref )


if __name__ == '__main__':

    test_tpsc_solver_h5()
