################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
# Author: H. U.R. Strand
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

from common import *

from triqs_tprf.linalg import inverse_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.utilities import G2_loc_fixed_fermionic_window_python
from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.bse import solve_lattice_bse

# -- Solve the lattice BSE for several fermionic window sizes
for nwf in [8, 10, 12, 20]:

    with HDFArchive('data_g2.h5', 'r') as a: p = a['p']
    p.nwf = nwf
    
    # -- DMFT impurity vertex
    p.chi_m = p.G2_iw_ph[('up','up')] - p.G2_iw_ph[('up','do')]
    p.chi_m = G2_loc_fixed_fermionic_window_python(p.chi_m, nwf=p.nwf)
    p.chi0_m = chi0_from_gg2_PH(p.G_w['up'], p.chi_m)
    p.gamma_m = inverse_PH(p.chi0_m) - inverse_PH(p.chi_m)

    # -- Lattice BSE
    g_wk = lattice_dyson_g_wk(mu=p.mu, e_k=p.e_k, sigma_w=p.sigma_w)[0:1, 0:1]
    p.chi_kw, p.chi0_kw = solve_lattice_bse(g_wk, p.gamma_m)

    with HDFArchive('data_bse_nwf{:03d}.h5'.format(nwf), 'w') as a: a['p'] = p
