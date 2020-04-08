"""
Read the local interaction tensor from UIJKL
and compare to the momentum dependent interaction UIJKL_Q after summing over Q.

Author: H. U.R. Strand
"""

import numpy as np
from vasp_crpa_parsers import read_vasp_crpa_momentum_space_interaction_to_triqs

from ase.calculators.vasp import Vasp2

calc = Vasp2(restart=True, directory='scf')

atoms = calc.get_atoms()

u_q = read_vasp_crpa_momentum_space_interaction_to_triqs('./crpa', atoms.cell, calc.kpts, verbose=False)

print(u_q)
