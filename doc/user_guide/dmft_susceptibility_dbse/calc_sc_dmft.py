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


p = ParameterCollection(
    U = 2.3,
    J_over_U = 0.4/2.3,
    N_target=4.0,
    n_k = 16,
    mu = +5.4,
    mu_min = +5.2,
    mu_max = +5.8,
    sc_iter_max = 10,
    mu_iter_max = 0,
    G_tol = 1e-4,
    N_tol = 5e-2,
    # DLR
    w_max = 4.,
    eps = 1e-8,
    )

p.solve = ParameterCollection(
    length_cycle = 100,
    n_warmup_cycles = int(1e6),
    n_cycles = int(1e9),
    )

p.init = ParameterCollection(
    beta = 25.0,
    n_iw = 400,
    n_tau = 4000,
    )

filename = f'./data/data_sc.h5' 
p0 = setup_dmft_calculation(p)

ps = solve_self_consistent_dmft(p0, verbose=False)
#ps = solve_self_consistent_dmft_fix_N(p0, verbose=False, filename=filename)

if mpi.is_master_node():
   with HDFArchive(filename, 'w') as a:
      a['ps'] = ps
