################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 by Hugo U. R. Strand
# Author: Hugo U. R. Strand
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

filename = './data/data_sc.h5'

if mpi.is_master_node():
    with HDFArchive(filename, 'r') as a:
        p = a['ps'][-1]
else: p = None
p = mpi.bcast(p)

p.solve.worm = True
p.solve.measure_G_l = False
p.solve.n_cycles = int(1e7)
p.solve.n_warmup_cycles = int(1e6)
p.solve.cfg_qmc = dict(
    WormEta=1, WormSearchEta=1, WormMeasGtau=1,
    PercentageWormInsert=0.3, PercentageWormReplace=0.1,
    Nwarmups2Plus=int(1e6),
    )

from w2dyn_cthyb import Solver
cthyb = Solver(**p.init.dict())
for bidx, g0 in cthyb.G0_iw: g0 << p.G0_w[bidx]
cthyb.solve(**p.solve.dict())
p.G_tau = cthyb.G_tau

if mpi.is_master_node():
    with HDFArchive(f'./data/data_g.h5', 'w') as a: a['p'] = p
