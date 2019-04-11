import os
import numpy as np

from ase.calculators.vasp import Vasp2

print '='*72
print '--> cRPA'
print '='*72

os.system('cp -r scf_exact crpa')

wannier90_win = """
num_wann =    3
num_bands=   96

# V(3d)-t2g states (bands 21-23)

exclude_bands = 1-20, 24-96

begin projections
 V:dxy;dxz;dyz
end projections
"""

with open('./crpa/wannier90.win', 'w') as fd:
    fd.write(wannier90_win)

calc = Vasp2(
    pp = 'PBE',
    setups = {'Sr' : '_sv_GW', 'V' : '_sv_GW', 'O' : '_GW'},
    command='mpirun vasp_std',
    restart = True,
    algo = 'CRPA',
    ncshmem = 1,
    precfock = 'Fast',
    ntarget_states = [1, 2, 3],
    lwrite_wanproj = True,
    lcrpaqout = True,
    loptics = None,
    nelm = None,
    #kpar = None,
    directory='crpa',
    txt = False)

try:
    calc.calculate(calc.atoms)
except:
    print '--> No energy from cRPA'
