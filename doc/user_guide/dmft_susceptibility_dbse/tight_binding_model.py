################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 by H. U.R. Strand
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

import numpy as np

from itertools import product

from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.wannier90 import parse_hopping_from_wannier90_hr_dat as read_hopping
from triqs_tprf.wannier90 import parse_lattice_vectors_from_wannier90_wout as read_units


def tight_binding_model(seed='sro', path='.', mu=7.832164838532275):

    units = read_units('sro.wout')
    hopping, num_wann = read_hopping('sro_hr.dat')

    hopping[(0,0,0)] -= np.eye(3) * mu
    hopping = { vec : np.kron(np.eye(2), mat) for vec, mat in hopping.items() } 
    
    orbital_names = [f'{s}_{o}' for s, o in product(['up', 'do'], range(num_wann))]
    
    tb_lattice = TBLattice(units=units, hopping=hopping,
        orbital_names=orbital_names, orbital_positions=[(0,0,0)]*2*num_wann)
    
    return tb_lattice
