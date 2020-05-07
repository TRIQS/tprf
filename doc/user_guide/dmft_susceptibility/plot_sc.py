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
from triqs.plot.mpl_interface import oplot, oplotr, oploti, plt

with HDFArchive('data_sc.h5', 'r') as a: ps = a['ps']
p = ps.objects[-1]
    
plt.figure(figsize=(3.25*2, 5))
subp = [2, 2, 1]

plt.subplot(*subp); subp[-1] += 1
plt.plot(ps.iter, ps.dG_l, 's-')
plt.ylabel('$\max | \Delta G_l |$')
plt.xlabel('Iteration')
plt.semilogy([], [])

plt.subplot(*subp); subp[-1] += 1
for b, g in p.G_l: p.G_l[b].data[:] = np.abs(g.data)
oplotr(p.G_l['up'], 'o-', label=None)
plt.ylabel('$| G_l |$')
plt.semilogy([], [])

plt.subplot(*subp); subp[-1] += 1
oplotr(p.G_tau_raw['up'], alpha=0.75, label='Binned')
oplotr(p.G_tau['up'], label='Legendre')
plt.legend(loc='best', fontsize=8)
plt.ylabel(r'$G(\tau)$')

plt.subplot(*subp); subp[-1] += 1
oploti(p.sigma_w[0,0], '.-', label=None)
plt.ylabel(r'$\Sigma(i\omega_n)$')
plt.xlim([-100, 100])

plt.tight_layout()
plt.savefig('figure_sc.svg')
plt.show()
