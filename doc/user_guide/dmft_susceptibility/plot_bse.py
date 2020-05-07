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
from triqs_tprf.lattice_utils import k_space_path

def get_path():    

    G = np.array([0.0, 0.0, 0.0]) * 2.*np.pi
    X = np.array([0.5, 0.5, 0.0]) * 2.*np.pi
    M = np.array([0.5, 0.0, 0.0]) * 2.*np.pi

    paths = [(G, X), (X, M), (M, G)]
    k_vecs, k_plot, K_plot = k_space_path(paths, num=100)
    K_labels = [r'$\Gamma$',r'$X$',r'$M$',r'$\Gamma$']

    return k_vecs, k_plot, K_plot, K_labels

def plot_chi_k(p, color=None):

    chi_k = p.chi_kw[:, Idx(0)]

    k_vecs, k_plot, K_plot, K_labels = get_path()
    kx, ky, kz = k_vecs.T

    interp = np.vectorize(lambda kx, ky, kz : np.squeeze(chi_k([kx, ky, kz]).real))
    interp = interp(kx, ky, kz)
    p.chi_G = np.squeeze(chi_k[Idx(0, 0, 0)].real)    

    line = plt.plot(k_plot, interp, '-', label=r'$N_{\nu} = $'+'${:d}$'.format(p.nwf))

    plt.gca().set_xticks(K_plot); plt.gca().set_xticklabels(K_labels)
    plt.grid(True); plt.ylabel(r'$\chi(\mathbf{Q})$')
    plt.legend(loc='best', fontsize=8)
    return line[0].get_color()

plt.figure(figsize=(3.25*2, 3))

ps, color = [], None
for filename in np.sort(glob.glob('data_bse_nwf*.h5')):
    with HDFArchive(filename, 'r') as a: p = a['p']
    subp = [1, 2, 1]
    plt.subplot(*subp); subp[-1] += 1
    color = plot_chi_k(p)
    plt.subplot(*subp); subp[-1] += 1
    plt.plot(1./p.nwf, p.chi_G, 'o', color=color)
    ps.append(p)

# -- Extrapolation to nwf -> oo
ps = ParameterCollections(objects=ps)
x, y = 1./ps.nwf, ps.chi_G
sidx = np.argsort(x)
x, y = x[sidx], y[sidx]
p = np.polyfit(x, y, 1)
y0 = np.polyval(p, 0)
X = np.linspace(0, x.max())
Y = np.polyval(p, X)

subp = [1, 2, 1]
plt.subplot(*subp); subp[-1] += 1
plt.plot(0, y0, 'rx')
plt.plot(0, 0.3479, 'r+')
plt.subplot(*subp); subp[-1] += 1
plt.plot(X, Y, '--k', lw=1.0, zorder=-100)
plt.plot(0, 0.3479, 'r+', label=r'Field')
plt.plot(0, y0, 'rx', label=r'BSE')
plt.grid(True); plt.legend(loc='best', fontsize=8)
plt.xlabel(r'$1/N_\nu$'); plt.ylabel(r'$\chi(\mathbf{0})$')
plt.title(
    r'$\lim_{n_\nu \rightarrow \infty} \, \chi(\mathbf{0}) \approx $' + \
    '${:3.4f}$'.format(y0))
plt.tight_layout()
plt.savefig('figure_bse.svg')
plt.show()
