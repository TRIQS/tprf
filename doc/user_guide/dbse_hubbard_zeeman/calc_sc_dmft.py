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

#from common import *

import triqs.utility.mpi as mpi
from common import setup_dmft_calculation, solve_self_consistent_dmft
from h5 import HDFArchive
from sys import argv

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections

p = ParameterCollection(
    t = 1.,
    t_prime = 0.,
    B = 0.,
    crystal_field = 0.,
    U = 0.5,
    N_target=1.0,
    n_l = 30,
    n_k = 8,
    mu = -0.8,
    sc_iter_max = 10,
    G_l_tol = 1e-4,
    )

p.solve = ParameterCollection(
    length_cycle = 15,
    n_warmup_cycles = int(1e7), #1e7
    n_cycles = int(1e8), #1e8
    )

p.init = ParameterCollection(
    beta = 5,
    n_iw = 1200,
    n_tau = 12000,
    )

filename = f'data_sc.h5'
p0 = setup_dmft_calculation(p)
ps = solve_self_consistent_dmft(p0, verbose=False)
if mpi.is_master_node():
   with HDFArchive(filename, 'w') as a:
      a['ps'] = ps
