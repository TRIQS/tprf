# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation for a Hubbard atom.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

from pytriqs.gf import Gf, Idx
from pytriqs.operators import c, c_dag
from h5 import HDFArchive
from pytriqs.utility import mpi # needed for pomerol2triqs

# ----------------------------------------------------------------------

from pomerol2triqs import PomerolED
from triqs_tprf.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
def make_calc():

    # ------------------------------------------------------------------
    # -- Hamiltonian

    p = ParameterCollection(
        beta = 0.5,
        U = 0.5,
        nw = 1,
        nwf = 15,
        V = 1.0,
        eps = 0.2,
        )

    p.nwf_gf = 4 * p.nwf
    p.mu = 0.5*p.U

    # ------------------------------------------------------------------

    ca_up, cc_up = c('0', 0), c_dag('0', 0)
    ca_do, cc_do = c('0', 1), c_dag('0', 1)

    ca0_up, cc0_up = c('1', 0), c_dag('1', 0)
    ca0_do, cc0_do = c('1', 1), c_dag('1', 1)

    docc = cc_up * ca_up * cc_do * ca_do
    nA = cc_up * ca_up + cc_do * ca_do
    hybridiz = p.V * (cc0_up * ca_up + cc_up * ca0_up + cc0_do * ca_do + cc_do * ca0_do)
    bath_lvl = p.eps * (cc0_up * ca0_up + cc0_do * ca0_do)

    p.H_int = p.U * docc
    p.H = -p.mu * nA + p.H_int + hybridiz + bath_lvl

    # ------------------------------------------------------------------
    # -- Exact diagonalization

    # Conversion from TRIQS to Pomerol notation for operator indices
    # TRIQS:   block_name, inner_index
    # Pomerol: site_label, orbital_index, spin_name
    index_converter = {
        ('0', 0) : ('loc', 0, 'up'),
        ('0', 1) : ('loc', 0, 'down'),
        ('1', 0) : ('loc', 1, 'up'),
        ('1', 1) : ('loc', 1, 'down'),
        }

    # -- Create Exact Diagonalization instance
    ed = PomerolED(index_converter, verbose=True)
    ed.diagonalize(p.H) # -- Diagonalize H

    p.gf_struct = [['0', [0, 1]]]

    # -- Single-particle Green's functions
    p.G_iw = ed.G_iw(p.gf_struct, p.beta, n_iw=p.nwf_gf)['0']

    # -- Particle-particle two-particle Matsubara frequency Green's function
    opt = dict(
        beta=p.beta, gf_struct=p.gf_struct,
        blocks=set([("0", "0")]),
        n_iw=p.nw, n_inu=p.nwf)

    p.G2_iw_ph = ed.G2_iw_inu_inup(channel='PH', **opt)[('0', '0')]

    filename = 'data_pomerol.h5'
    with HDFArchive(filename,'w') as res:
        res['p'] = p

    import os
    os.system('tar czvf data_pomerol.tar.gz data_pomerol.h5')
    os.remove('data_pomerol.h5')
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc()
