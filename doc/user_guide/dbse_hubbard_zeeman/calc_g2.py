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

import itertools

import numpy as np

import triqs.utility.mpi as mpi
from h5 import HDFArchive

from triqs.gfs import Gf, MeshImFreq, Fourier, LegendreToMatsubara, BlockGf, inverse, Idx, MeshProduct
from triqs.gfs.tools import fit_legendre

from triqs.operators import n, c, c_dag, Operator
from triqs.operators.util.hamiltonians import h_int_kanamori, h_int_density
from triqs.operators.util.op_struct import set_operator_structure

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections

import sys

if mpi.is_master_node():
    with HDFArchive(f'data_sc.h5', 'r') as a:
        p = a['ps'][-1]
else: p = None
p = mpi.bcast(p)

# -- Rewrite problem with completely spin-orbital diagonal Green's function
# -- and relabel local interaction operators.

num_orbitals = 1
spin_names = ('up','do')
orb_names = list(range(num_orbitals))
p.init.gf_struct = set_operator_structure(spin_names, orb_names, off_diag=False)

from pyed.OperatorUtils import relabel_operators

from_list = [ op(spin, orb) for spin, orb, op in itertools.product(spin_names, orb_names, [c, c_dag]) ]
to_list = [ op(spin+f'_{orb}', 0)  for spin, orb, op in itertools.product(spin_names, orb_names, [c, c_dag]) ]

p.solve.h_int = relabel_operators(p.solve.h_int, from_list, to_list)

WormComponents = [ list(x) for x in itertools.product(range(2*num_orbitals), repeat=4) ] 

print(f'WormComponents =\n{WormComponents}')

cfg_qmc=dict(
    PercentageWormInsert=0.3,
    PercentageWormReplace=0.1,
    WormEta=1,
    WormSearchEta=1,
    Nwarmups2Plus=10000,
    WormMeasG4iw=1,
    FourPnt=8,
    WormPHConvention=0, # Important to get correct frequency structure
    N4iwb=30,
    N4iwf=30,
    WormComponents=WormComponents, # Add G2 components to sample here.
    )

p.solve.n_warmup_cycles = int(1e7) 
p.solve.n_cycles = int(1e8) 
p.solve.measure_G_l = False

p.solve.worm = True
p.solve.cfg_qmc = cfg_qmc

from w2dyn_cthyb import Solver

cthyb = Solver(**p.init.dict())

cthyb.G0_iw['up_0'][0, 0] << p.G0_w['up'][0, 0]
cthyb.G0_iw['do_0'][0, 0] << p.G0_w['do'][0, 0]

cthyb.solve(**p.solve.dict())
p.G2_worm_components = cthyb.G2_worm_components

if mpi.is_master_node():
    with HDFArchive(f'data_g2.h5', 'w') as a: a['p'] = p
