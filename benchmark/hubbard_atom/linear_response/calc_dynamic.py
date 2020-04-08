# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation 
    for a Hubbard atom with two bath sites. 

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

"""

# ----------------------------------------------------------------------

import os
import numpy as np

from pytriqs.operators import c, c_dag
from h5 import HDFArchive
from pytriqs.gf import GfImTime, GfImFreq, BlockGf

# ----------------------------------------------------------------------

from pytriqs.utility import mpi

# ----------------------------------------------------------------------

from triqs_tprf.analytic_hubbard_atom import analytic_hubbard_atom
from triqs_tprf.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
def make_calc(beta=2.0, nwf=8):
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    p = ParameterCollection(
        beta = beta,
        U = 5.0,
        nw = 1,
        nwf = nwf,
        nwf_gf = 2*nwf,
        )

    ana = analytic_hubbard_atom(**p.dict())

    p.chi = np.sum(ana.chi_m.data) / p.beta**2
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_dynamic_beta%6.6f_nwf%i.h5' % (p.beta, p.nwf)
    with HDFArchive(filename,'w') as res:
        res['p'] = p
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    beta_vec = np.array([0.01, 0.1, 1.0, 10., 100., 200., 500., 1000.])

    for beta in beta_vec:
        path = 'dynamic_beta%6.6f' % beta
        print('--> path:', path)
        os.mkdir(path)
        os.chdir(path)

        for nwf in [8, 16, 32, 64, 128, 256]:
            make_calc(beta=beta, nwf=nwf)

        os.chdir('../')
