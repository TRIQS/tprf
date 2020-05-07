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

ps = [] 
filenames = np.sort(glob.glob('data_B_*.h5'))
for filename in filenames:
    with HDFArchive(filename, 'r') as a:
        p = a['ps'].objects[-1]
    ps.append(p)
ps = ParameterCollections(ps)

B, M = ps.B, ps.M
B = np.concatenate((-B[1:][::-1], B))
M = np.concatenate((-M[1:][::-1], M))
p = np.polyfit(M, B, 5)
m = np.linspace(-0.5, 0.5, num=1000)
b = np.polyval(p, m)
chi = 1./np.polyval(np.polyder(p, 1), 0.).real

plt.figure(figsize=(3.25*1.5, 2.5*1.5))
plt.title(r'$\chi = \frac{dM}{dB}|_{B=0} \approx $'+'$ {:3.4f}$'.format(chi))
plt.plot(B, M, 'o', alpha=0.75, label='DMFT field')
plt.plot(b, m, '-', alpha=0.75, label='Poly fit')
plt.legend(loc='upper left'); plt.grid(True)
plt.xlabel(r'$B$'); plt.ylabel(r'$M$')
plt.xlim([B.min(), B.max()]); plt.ylim([-0.5, 0.5])
plt.tight_layout()
plt.savefig('figure_field.svg')
plt.show()
