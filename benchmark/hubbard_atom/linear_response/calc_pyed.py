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

from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
def make_calc(beta=2.0, h_field=0.0):
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    p = ParameterCollection(
        beta = beta,
        h_field = h_field,
        U = 5.0,
        ntau = 40,
        niw = 15,
        )

    p.mu = 0.5*p.U
    
    # ------------------------------------------------------------------

    print '--> Solving SIAM with parameters'
    print p
    
    # ------------------------------------------------------------------

    up, do = 'up', 'dn'
    docc = c_dag(up,0) * c(up,0) * c_dag(do,0) * c(do,0)
    mA = c_dag(up,0) * c(up,0) - c_dag(do,0) * c(do,0)
    nA = c_dag(up,0) * c(up,0) + c_dag(do,0) * c(do,0)

    p.H = -p.mu * nA + p.U * docc + p.h_field * mA
    
    # ------------------------------------------------------------------

    fundamental_operators = [c(up,0), c(do,0)]
    
    ed = TriqsExactDiagonalization(p.H, fundamental_operators, p.beta)

    g_tau = GfImTime(beta=beta, statistic='Fermion', n_points=40, indices=[0])
    g_iw = GfImFreq(beta=beta, statistic='Fermion', n_points=10, indices=[0])

    p.G_tau = BlockGf(name_list=[up,do], block_list=[g_tau]*2, make_copies=True)
    p.G_iw = BlockGf(name_list=[up,do], block_list=[g_iw]*2, make_copies=True)
    
    ed.set_g2_tau(p.G_tau[up], c(up,0), c_dag(up,0))
    ed.set_g2_tau(p.G_tau[do], c(do,0), c_dag(do,0))

    ed.set_g2_iwn(p.G_iw[up], c(up,0), c_dag(up,0))
    ed.set_g2_iwn(p.G_iw[do], c(do,0), c_dag(do,0))

    p.magnetization = ed.get_expectation_value(0.5 * mA)
    p.magnetization2 = ed.get_expectation_value(0.25 * mA * mA)
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_pyed_h_field_%4.4f.h5' % h_field
    with HDFArchive(filename,'w') as res:
        res['p'] = p
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    beta_vec = np.array([0.01, 0.1, 1.0, 10., 100., 200., 500., 1000.])
    
    for beta in beta_vec:
        path = 'pyed_beta%6.6f' % beta
        print '-->path: ', path
        os.mkdir(path)
        os.chdir(path)

        #h_field_vec = np.linspace(-10., 10., num=21) / np.sqrt(beta)
        h_field_vec = np.linspace(-1., 1., num=21) / beta
        
        for h_field in h_field_vec:
            make_calc(beta=beta, h_field=h_field)

        os.chdir('../')
