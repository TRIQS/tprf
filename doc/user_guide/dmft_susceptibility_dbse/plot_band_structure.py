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

from tight_binding_model import tight_binding_model
t_r = tight_binding_model()

G, M, X = [0., 0., 0.], [-0.25, 0.25, 0.25], [0., 0., 0.5]
labels = [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$', ]
paths = [(G, X), (X, M), (M, G),]

from triqs_tprf.lattice_utils import k_space_path
k_vecs, k_plot, K_plot = k_space_path(paths, bz=t_r.bz, num=32)

import matplotlib.pyplot as plt
plt.figure(figsize=(3.25*2, 3))

from numpy.linalg import eigvalsh as eigv
e_k_interp = [ eigv(t_r.tb.fourier(k)) for k in k_vecs ]
plt.plot(k_plot, e_k_interp, '-k')

e_k = t_r.fourier(t_r.get_kmesh(n_k=[16]*3))
e_k_interp_ref = [ eigv(e_k(tuple(k))) for k in k_vecs @ t_r.bz.units ]
plt.plot(k_plot, e_k_interp_ref, '-r', lw=0.5)

plt.xticks(ticks=K_plot, labels=labels)
plt.ylabel(r'$\epsilon(\mathbf{k})$')
plt.grid(True); plt.tight_layout();
plt.savefig('figure_sro_band_structure.svg')
plt.show()

