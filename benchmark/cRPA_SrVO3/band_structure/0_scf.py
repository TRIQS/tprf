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
    cell=np.diag([3.843]*3),
    )

#from ase.visualize import view
#view(atoms)

from ase.calculators.vasp import Vasp2

calc = Vasp2(
    system='SrVO3',
    command='mpirun vasp_std',
    pp = 'PBE',
    setups = {'Sr' : '_sv_GW', 'V' : '_sv_GW', 'O' : '_GW'},    
    gamma = True,
    ismear = 0,
    ediff = 1e-8,
    nbands = 36,
    kpts=(4, 4, 4), directory='scf', txt=False)

atoms.set_calculator(calc)
atoms.get_potential_energy()
