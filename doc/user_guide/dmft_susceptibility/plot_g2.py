################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
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

from common import *
from triqs.plot.mpl_interface import oplot, oplotr, plt

with HDFArchive('data_g2.h5', 'r') as a: p = a['p']

plt.figure(figsize=(3.25*2.2, 2))
subp = [1, 3, 1]
opt = dict(cmap=plt.get_cmap('terrain_r'), vmin=0.0, vmax=0.01)
plt.subplot(*subp); subp[-1] += 1
plt.title(r'$\chi^{(0)}_m$')
plt.imshow(np.squeeze(p.chi0_m.data).real, **opt)
plt.colorbar()
plt.subplot(*subp); subp[-1] += 1
plt.title(r'$\chi_m$')
plt.imshow(np.squeeze(p.chi_m.data).real, **opt)
plt.colorbar()
plt.subplot(*subp); subp[-1] += 1
plt.title(r'$\Gamma_m - U$')
plt.imshow((np.squeeze(p.gamma_m.data) - p.U).real,
           cmap=plt.get_cmap('RdBu_r'), vmin=-5, vmax=5)
plt.colorbar()

plt.tight_layout()
plt.savefig('figure_g2.svg')
plt.show()
