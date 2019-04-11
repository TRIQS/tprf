from ase.build import bulk
from ase.calculators.vasp import Vasp2

import os

print '='*72, '\n', '--> Band structure\n', '='*72

os.system('cp -r scf bs')
calc = Vasp2(
    system='SrVO3',
    command='mpirun vasp_std',
    pp = 'PBE',
    setups = {'Sr' : '_sv_GW', 'V' : '_sv_GW', 'O' : '_GW'},
    isym = 0, icharg = 11, kpts = dict(path='GXMGRX', npoints=200),
    restart=True, directory='bs', txt=False)

calc.get_potential_energy()
 
