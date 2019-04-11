from ase.calculators.vasp import Vasp2

calc_load = Vasp2(restart=True, directory='bs')

bs = calc_load.band_structure() # ASE Band structure object
bs.plot(emin=-5)               # Plot the band structure
