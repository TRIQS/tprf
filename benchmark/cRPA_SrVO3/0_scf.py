import numpy as np
from ase import Atoms

atoms = Atoms(
    symbols=['Sr'] + ['V'] + ['O']*3,
    scaled_positions=[
        [0., 0., 0.], # Sr
        [0.5, 0.5, 0.5], # V 
        [0.5, 0.5, 0.], # O1
        [0., 0.5, 0.5], # O2
        [0.5, 0., 0.5], # O3
        ],
    cell=np.diag([3.84652]*3),
    )

from ase.calculators.vasp import Vasp2

print('='*72)
print('--> SCF')
print('='*72)

calc = Vasp2(
    command='mpirun vasp_std',
    system='SrVO3',
    pp = 'PBE',
    setups = {'Sr' : '_sv_GW', 'V' : '_sv_GW', 'O' : '_GW'},
    nbands = 36,
    gamma = True,
    ismear = 0,
    ediff = 1e-8,
    kpts=(2, 2, 2),
    kpar=2,
    directory='scf',
    txt = False)

atoms.set_calculator(calc)
atoms.get_potential_energy()
