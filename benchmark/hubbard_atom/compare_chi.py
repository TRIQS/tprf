
""" Compare analytic solution for chi_m with numerical ED results 

Author: Hugo U.R. Strand (2017), hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from triqs_tprf.hubbard_atom import chi_ph_magnetic
from triqs_tprf.hubbard_atom import gamma_ph_magnetic
from triqs_tprf.hubbard_atom import single_particle_greens_function

# ----------------------------------------------------------------------

from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
    
# ----------------------------------------------------------------------
class Dummy():
     def __init__(self):
        pass

# ----------------------------------------------------------------------
def read_pomerol_data(filename):

    data = Dummy()

    #print '--> Loading:', filename
    with HDFArchive(filename, 'r') as h5:
        for key, value in h5.items():
            setattr(data, key, value)                

    a, b =  data.G2_iw_ph[('up', 'up')].data.shape[:2]
    nw, nwf = (a+1)/2, b/2
    data.nw, data.nwf = nw, nwf    
            
    data.beta = data.params['beta']
    data.U = data.params['U']

    # -- Single particle Green's function
            
    g_mat = block_iw_AB_to_matrix_valued(data.G_iw)
    g_mat.name = 'g_mat'

    data.G_iw = data.G_iw['up']

    # -- Generalized susceptibility in magnetic PH channel

    data.chi_m = data.G2_iw_ph[('up', 'up')] - data.G2_iw_ph[('up','dn')]
    data.chi0_m = chi0_from_gg2_PH(g_mat, data.chi_m)

    data.label = r'Pomerol'

    return data
    
# ----------------------------------------------------------------------
def analytic_solution(beta, U, nw, nwf):

    data = Dummy()

    g_iw = single_particle_greens_function(beta=beta, U=U, nw=nw*2)
    data.G_iw = g_iw

    # make block gf of the single gf
    G_iw_block = BlockGf(name_list=['up', 'dn'], block_list=[g_iw, g_iw])
    g_mat = block_iw_AB_to_matrix_valued(G_iw_block)
    
    data.chi_m = chi_ph_magnetic(beta=beta, U=U, nw=nw, nwf=nwf)
    data.chi0_m = chi0_from_gg2_PH(g_mat, data.chi_m)

    data.label = r'Analytic'
    return data

# ----------------------------------------------------------------------
if __name__ == '__main__':

    pom = read_pomerol_data('data_pomerol.h5')
    ana = analytic_solution(beta=pom.beta, U=pom.U, nw=pom.nw, nwf=pom.nwf)

    # -- Compare single-particle Green's functions
    np.testing.assert_array_almost_equal(pom.G_iw.data, ana.G_iw.data)

    # -- Compare generalized susceptibility (magnetic channel)
    np.testing.assert_array_almost_equal(pom.chi_m.data, ana.chi_m.data)
    np.testing.assert_array_almost_equal(pom.chi0_m.data, ana.chi0_m.data)

    print 'ok! Analytic chi_m agrees with ED (pomerol).'
