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

if mpi.is_master_node():
    with HDFArchive('data_B_0.000000.h5', 'r') as a:
        p = a['ps'].objects[-1]
else: p = None
p = mpi.bcast(p)

# -- Sample G2

p.solve.n_cycles = int(1e9 / 40.)
p.solve.measure_G_l = False
p.solve.measure_G_tau = False
p.solve.measure_G2_iw_ph = True
p.solve.measure_G2_blocks = set([('up','up'), ('up','do')])
p.solve.measure_G2_n_bosonic = 1
p.solve.measure_G2_n_fermionic = 20

cthyb = triqs_cthyb.Solver(**p.init.dict())
cthyb.G0_iw << p.G0_w
cthyb.solve(**p.solve.dict())
p.G2_iw_ph = cthyb.G2_iw_ph.copy()

# -- Compute DMFT impurity vertex

from triqs_tprf.linalg import inverse_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH

p.chi_m = p.G2_iw_ph[('up','up')] - p.G2_iw_ph[('up','do')]
p.chi0_m = chi0_from_gg2_PH(p.G_w['up'], p.chi_m)
p.gamma_m = inverse_PH(p.chi0_m) - inverse_PH(p.chi_m)

del p.solve.measure_G2_blocks
if mpi.is_master_node():
    with HDFArchive('data_g2.h5', 'w') as a: a['p'] = p
