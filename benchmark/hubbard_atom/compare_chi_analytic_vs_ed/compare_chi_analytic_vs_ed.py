
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
def plot_res(ana, pom):

    from triqs_tprf.plot import plot_g2, get_g2_opt
    from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt

    opt = get_g2_opt(pom.chi, cplx='re', cut=0.1)
    
    plt.figure(figsize=(6, 6))
    plot_g2(pom.G2_iw_ph, cplx='re', opt=opt, idx_labels=[r'\uparrow', r'\downarrow'])

    plt.figure(figsize=(6, 6))
    plot_g2(pom.chi0, cplx='re', opt=opt, idx_labels=[r'\uparrow', r'\downarrow'])

    plt.figure(figsize=(6, 6))
    plot_g2(pom.chi, cplx='re', opt=opt, idx_labels=[r'\uparrow', r'\downarrow'])

    plt.figure(figsize=(3, 3))
    plot_g2(pom.chi_m, cplx='re', opt=opt, idx_labels=[r'\uparrow', r'\downarrow'])
    
    plt.figure(figsize=(3, 3))
    plot_g2(pom.chi0_m, cplx='re', opt=opt, idx_labels=[r'\uparrow', r'\downarrow'])

    opt = get_g2_opt(pom.gamma, cplx='re', cut=0.1)
    opt['vmin'] = -5.
    opt['vmax'] = +5.

    plt.figure(figsize=(6, 6))
    plot_g2(pom.gamma, cplx='re', opt=opt, idx_labels=[r'\uparrow', r'\downarrow'])
    
    plt.figure(figsize=(3, 3))
    plot_g2(ana.gamma_m, cplx='re', opt=opt, idx_labels=[r'\uparrow', r'\downarrow'])

    plt.show()
    
# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_pomerol.h5'
    with HDFArchive(filename, 'r') as h5:
         pom = h5['p']         
     
    ana = analytic_hubbard_atom(
         beta=pom.beta, U=pom.U,
         nw=pom.nw, nwf=pom.nwf, nwf_gf=pom.nwf_gf)

    # -- Compare single-particle Green's functions
    np.testing.assert_array_almost_equal(pom.G_iw[0,0].data, ana.G_iw[0,0].data)

    # -- Compare generalized susceptibility (magnetic channel)
    np.testing.assert_array_almost_equal(pom.chi_m.data, ana.chi_m.data)
    np.testing.assert_array_almost_equal(pom.chi0_m.data, ana.chi0_m.data)

    print 'ok! Analytic chi_m agrees with ED (pomerol).'

    plot_res(ana, pom)
    
    
