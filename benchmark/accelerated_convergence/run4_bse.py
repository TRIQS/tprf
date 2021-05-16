from common import *

from triqs_tprf.linalg import inverse_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.utilities import G2_loc_fixed_fermionic_window_python
from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.bse import solve_lattice_bse

from triqs.lattice import BravaisLattice,BrillouinZone
import triqs.gf as gf

# Single-point kmesh
BL = BravaisLattice(units = [(1,0,0) ,(0,1,0),(0,0,1) ]) #linear lattice
onepoint = gf.MeshBrillouinZone(BrillouinZone(BL), 1)

# -- Solve the lattice BSE for several fermionic window sizes
for nwf in [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20]:

    with HDFArchive('data_g2.h5', 'r') as a: p = a['p']
    p.nwf, p.tail_corr_nwf = nwf, 100
    
    # -- DMFT impurity vertex
    p.chi_m = p.G2_iw_ph[('up','up')] - p.G2_iw_ph[('up','do')]
    p.chi_m = G2_loc_fixed_fermionic_window_python(p.chi_m, nwf=p.nwf)
    p.chi0_m = chi0_from_gg2_PH(p.G_w['up'], p.chi_m)
    p.gamma_m = inverse_PH(p.chi0_m) - inverse_PH(p.chi_m)

    # -- Lattice BSE
    g_wk = lattice_dyson_g_wk(mu=p.mu, e_k=p.e_k, sigma_w=p.sigma_w)[0:1, 0:1]
    p.chi_kw, p.chi0_kw = solve_lattice_bse(g_wk, p.gamma_m, tail_corr_nwf=p.tail_corr_nwf)


    # We now calculate the same thing using the impurity Green's function
    # To use the usual lattice BSE code, we need to promote p.g_w to a g_wk object
    # With a single k point
    aux_mesh = gf.MeshProduct(p.g_w.mesh,onepoint)
    g_imp = gf.Gf(indices=[0,],mesh=aux_mesh)
    
    for w,k in g_imp.mesh:
        g_imp[w,k][0,0] = p.g_w[w][0,0]
    
    p.chi_w_imp, p.chi0_w_imp = solve_lattice_bse(g_imp, p.gamma_m, tail_corr_nwf=p.tail_corr_nwf)

    with HDFArchive('data_bse_nwf{:03d}.h5'.format(nwf), 'w') as a: a['p'] = p
