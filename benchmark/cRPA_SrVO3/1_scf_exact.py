
import os
import numpy as np

from ase.calculators.vasp import Vasp2

print '='*72
print '--> SCF exact'
print '='*72

os.system('cp -r scf scf_exact')

calc = Vasp2(
    command='mpirun vasp_std',
    pp = 'PBE',
    setups = {'Sr' : '_sv_GW', 'V' : '_sv_GW', 'O' : '_GW'},
    restart = True,
    algo = 'Exact',
    nelm = 1,
    loptics = True,
    ediff = None,
    nbands = 112,
    directory='scf_exact',
    txt = False)

calc.get_potential_energy()
