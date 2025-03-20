################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 by Hugo U. R. Strand
# Author: Hugo U. R. Strand
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

from triqs_tprf.wannier90 import parse_hopping_from_wannier90_hr_dat
from triqs_tprf.wannier90 import parse_lattice_vectors_from_wannier90_wout


def tight_binding_model(mu=11.6591, path='./calc_dft/wannier_fit/'):

    units = parse_lattice_vectors_from_wannier90_wout(path + 'sro.wout')
    hopping, n_orb = parse_hopping_from_wannier90_hr_dat(path + 'sro_hr.dat')

    hopping[(0,0,0)] -= np.eye(n_orb) * mu
    hopping = { vec : np.kron(np.eye(2), mat) for vec, mat in hopping.items() } 
    orbital_names = [f'{s}_{o}' for s, o in product(['up', 'do'], range(n_orb))]
    
    tb_lattice = TBLattice(units=units, hopping=hopping,
        orbital_names=orbital_names, orbital_positions=[(0,0,0)]*2*n_orb)
    
    return tb_lattice


def dft_band_structure(mu=11.6591, filename='./calc_dft/band_structure/sro_bands.dat.gnu'):

    data = np.loadtxt(filename)
    
    n_k = len(np.unique(data[:, 0]))
    k = data[:n_k, 0]
    e_k = data[:, 1].reshape(-1, n_k)

    K_plot = [0.0000, 0.8961, 1.7923, 3.0596]
    labels = [r'$\Gamma$', '$X$', '$M$', r'$\Gamma$']

    return k, e_k.T - mu, K_plot


if __name__ == '__main__':

    t_r = tight_binding_model()
    print(t_r.orbital_names)

    G, M, X = [0., 0., 0.], [0., 0., 0.5], [-0.25, 0.25, 0.25]
    labels = [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$', ]
    paths = [(G, X), (X, M), (M, G),]

    from triqs.lattice.utils import k_space_path
    k_vecs, k_plot, K_plot = k_space_path(paths, bz=t_r.bz, num=128, relative_coordinates=True)

    from numpy.linalg import eigvalsh as eigv
    e_k_interp = [ eigv(t_r.tb.fourier(k)) for k in k_vecs ]

    k, e_k_qe, K_plot_qe = dft_band_structure()
    k *= k_plot[-1] / k[-1]
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(3.25*2, 6))

    plt.plot(k, e_k_qe, '-g', alpha=0.5)
    plt.plot(k_plot, e_k_interp, '--b', alpha=0.5)
    
    plt.plot([], [], '-g', alpha=0.5, label='DFT QE@PAW+PBE')
    plt.plot([], [], '--b', alpha=0.5, label=r'Wannier Ru $4d$, $t_{2g}$')

    plt.xticks(ticks=K_plot, labels=labels)

    plt.ylim([-9, 7])
    plt.xlim([k_plot.min(), k_plot.max()])
    plt.ylabel(r'$\epsilon(\mathbf{k})$')
    plt.legend(loc='best')
    plt.grid(True); plt.tight_layout();
    plt.savefig('figure_sro_band_structure.svg')
    plt.show()

