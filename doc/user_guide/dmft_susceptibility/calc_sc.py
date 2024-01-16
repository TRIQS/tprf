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
    B = 0.5,
    U = 5.,
    #U = 1e-10,
    #J_over_U = 0.1,
    num_orbitals = 1,
    #mu = 2. + 5.,
    mu = 2. + 0.,
    n_k = 16,
    n_iter = 10,
    #n_iter = 1,
    G_tol = 1e-4,
    eps = 1e-10,
    #w_max = 20.0,
    w_max = 40.0,
    #w_max = 80.0,
    #w_max = 200.0,
    )

p.solve = ParameterCollection(
    length_cycle = 20,
    #n_warmup_cycles = int(1e4),
    #n_warmup_cycles = int(1e5),
    n_warmup_cycles = int(1e6),
    #n_cycles = int(8e6),
    #n_cycles = int(1e5),
    n_cycles = int(1e6),
    #n_cycles = int(1e7),
    #n_cycles = int(1e8),
    move_double = False,
    measure_G_l = False,
    use_norm_as_weight = True,
    measure_density_matrix = True,
    perform_post_proc = True,
    perform_tail_fit = True,    
    )

p.init = ParameterCollection(
    beta = 1. + 1/3,
    n_iw = 400,
    n_tau = 1024*16,
    )

filename = './data/data_sc_start.h5'

import os
if os.path.isfile(filename):
    if mpi.is_master_node(): print(f'--> Loading: {filename}')
    with HDFArchive(filename, 'r') as a:
        p0 = a['ps'][-1]

    p0.G_tol = 1e-5
    p0.solve.n_cycles = int(1e7)
    
else:
    p0 = setup_dmft_calculation(p)
    
ps = solve_self_consistent_dmft(p0)

if mpi.is_master_node():
    with HDFArchive('./data/data_sc.h5', 'w') as a:
        a['ps'] = ParameterCollections(ps)
