
""" Compare analytic solution for chi_m with numerical ED results 

Author: Hugo U.R. Strand (2017), hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.analytic_hubbard_atom import analytic_hubbard_atom

# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_pomerol.h5'
    with HDFArchive(filename, 'r') as h5:
         pom = h5['p']         
     
    ana = analytic_hubbard_atom(
         beta=pom.beta, U=pom.U,
         nw=pom.nw, nwf=pom.nwf, nwf_gf=pom.nwf_gf)
    
    # -- Compare single-particle Green's functions
    np.testing.assert_array_almost_equal(pom.G_iw.data, ana.G_iw.data)

    # -- Compare generalized susceptibility (magnetic channel)
    np.testing.assert_array_almost_equal(pom.chi_m.data, ana.chi_m.data)
    np.testing.assert_array_almost_equal(pom.chi0_m.data, ana.chi0_m.data)

    print 'ok! Analytic chi_m agrees with ED (pomerol).'
