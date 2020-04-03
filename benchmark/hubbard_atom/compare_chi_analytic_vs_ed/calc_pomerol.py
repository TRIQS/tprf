# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation for a Hubbard atom.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

from pytriqs.gf import Gf
from pytriqs.operators import c, c_dag
from h5 import HDFArchive
from pytriqs.utility import mpi # needed for pomerol2triqs

# ----------------------------------------------------------------------

from pomerol2triqs import PomerolED

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH
from triqs_tprf.linalg import inverse_PH

# ----------------------------------------------------------------------
def make_calc():
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    p = ParameterCollection(
        beta = 1.0,
        U = 5.0,
        nw = 1,
        nwf = 20,
        )

    p.nwf_gf = 4 * p.nwf
    p.mu = 0.5*p.U
    
    # ------------------------------------------------------------------

    ca_up, cc_up = c('0', 0), c_dag('0', 0)
    ca_do, cc_do = c('0', 1), c_dag('0', 1)
    
    docc = cc_up * ca_up * cc_do * ca_do
    nA = cc_up * ca_up + cc_do * ca_do
    p.H = -p.mu * nA + p.U * docc
    
    # ------------------------------------------------------------------
    # -- Exact diagonalization

    # Conversion from TRIQS to Pomerol notation for operator indices
    # TRIQS:   block_name, inner_index
    # Pomerol: site_label, orbital_index, spin_name
    index_converter = {
        ('0', 0) : ('loc', 0, 'up'),
        ('0', 1) : ('loc', 0, 'down'),
        }

    # -- Create Exact Diagonalization instance
    ed = PomerolED(index_converter, verbose=True)
    ed.diagonalize(p.H) # -- Diagonalize H

    gf_struct = [['0', [0, 1]]]

    # -- Single-particle Green's functions
    p.G_iw = ed.G_iw(gf_struct, p.beta, n_iw=p.nwf_gf)['0']

    # -- Particle-particle two-particle Matsubara frequency Green's function
    opt = dict(
        beta=p.beta, gf_struct=gf_struct,
        blocks=set([("0", "0")]),
        n_iw=p.nw, n_inu=p.nwf)
    
    p.G2_iw_ph = ed.G2_iw_inu_inup(channel='PH', **opt)[('0', '0')]
    
    # ------------------------------------------------------------------
    # -- Generalized susceptibility in magnetic PH channel
            
    p.chi_m = Gf(mesh=p.G2_iw_ph.mesh, target_shape=[1, 1, 1, 1])    
    p.chi_m[0, 0, 0, 0] = p.G2_iw_ph[0, 0, 0, 0] - p.G2_iw_ph[0, 0, 1, 1]
    
    p.chi0_m = chi0_from_gg2_PH(p.G_iw, p.chi_m)
    p.label = r'Pomerol'
    
    # ------------------------------------------------------------------
    # -- Generalized susceptibility in PH channel
            
    p.chi = chi_from_gg2_PH(p.G_iw, p.G2_iw_ph)
    p.chi0 = chi0_from_gg2_PH(p.G_iw, p.G2_iw_ph)
    p.gamma = inverse_PH(p.chi0) - inverse_PH(p.chi)

    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_pomerol.h5'
    with HDFArchive(filename,'w') as res:
        res['p'] = p
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc()
