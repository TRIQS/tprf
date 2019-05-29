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

p = ParameterCollection(
    t = 1.,
    B = 0.,
    U = 10.,
    mu = 0.,
    n_k = 16,
    n_iter = 10,
    G_l_tol = 1e-5,
    )

p.solve = ParameterCollection(
    length_cycle = 10,
    n_warmup_cycles = 10000,
    n_cycles = int(8e6),
    move_double = False,
    measure_G_l = True,
    )

p.init = ParameterCollection(
    beta = 1.,
    n_l = 10,
    n_iw = 400,
    n_tau = 4000,
    gf_struct = [('up',[0]), ('do',[0])])

p0 = setup_dmft_calculation(p)
ps = solve_self_consistent_dmft(p0)

if mpi.is_master_node():
    with HDFArchive('data_sc.h5', 'w') as a:
        a['ps'] = ParameterCollections(ps)
