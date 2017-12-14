# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation for a Hubbard atom.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

from pytriqs.operators import c, c_dag
from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi # needed for pomerol2triqs

# ----------------------------------------------------------------------

from pomerol2triqs import PomerolED

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH

# ----------------------------------------------------------------------
def make_calc():
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    p = ParameterCollection(
        beta = 2.0,
        U = 5.0,
        nw = 2,
        nwf = 3,
        )

    p.nwf_gf = 4
    p.mu = 0.5*p.U
    
    # ------------------------------------------------------------------

    up, do = 'up', 'dn'
    docc = c_dag(up,0) * c(up,0) * c_dag(do,0) * c(do,0)
    nA = c_dag(up,0) * c(up,0) + c_dag(do,0) * c(do,0)
    p.H = -p.mu * nA + p.U * docc
    
    # ------------------------------------------------------------------
    # -- Exact diagonalization

    # Conversion from TRIQS to Pomerol notation for operator indices
    # TRIQS:   block_name, inner_index
    # Pomerol: site_label, orbital_index, spin_name
    index_converter = {
        (up, 0) : ('loc', 0, 'up'),
        (do, 0) : ('loc', 0, 'down'),
        }

    # -- Create Exact Diagonalization instance
    ed = PomerolED(index_converter, verbose=True)
    ed.diagonalize(p.H) # -- Diagonalize H

    gf_struct = {up : [0], do : [0]}

    # -- Single-particle Green's functions
    p.G_iw = ed.G_iw(gf_struct, p.beta, n_iw=p.nwf_gf)

    # -- Particle-particle two-particle Matsubara frequency Green's function
    opt = dict(
        beta=p.beta, gf_struct=gf_struct,
        blocks=set([("up", "up"), ("up", "dn")]),
        n_iw=p.nw, n_inu=p.nwf)
    
    G2_iw_ph = ed.G2_iw_inu_inup(channel='PH', **opt)
    
    # ------------------------------------------------------------------
    # -- Generalized susceptibility in magnetic PH channel
            
    g_mat = block_iw_AB_to_matrix_valued(p.G_iw)
    g_mat.name = 'g_mat'
    p.G_iw = p.G_iw['up'] # only store one component of G_iw

    p.chi_m = G2_iw_ph[('up', 'up')] - G2_iw_ph[('up','dn')]
    p.chi0_m = chi0_from_gg2_PH(g_mat, p.chi_m)
    p.label = r'Pomerol'
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_pomerol.h5'
    with HDFArchive(filename,'w') as res:
        res['p'] = p
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc()
