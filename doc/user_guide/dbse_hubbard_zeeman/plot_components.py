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
    M = np.array([0.5, 0.5, 0.0]) * 2.*np.pi
    X = np.array([0.5, 0.0, 0.0]) * 2.*np.pi

    paths = [(G, M), (M, X), (X, G)]
    k_vecs, k_plot, K_plot = k_space_path(paths, num=100)
    K_labels = [r'$\Gamma$',r'$M$',r'$X$',r'$\Gamma$']

    return k_vecs, k_plot, K_plot, K_labels

def plot_chi_k(p, chi_kw, w, component, color=None):

    chi_k = chi_kw[:, Idx(w)][component]

    k_vecs, k_plot, K_plot, K_labels = get_path()
    kx, ky, kz = k_vecs.T

    interp = np.vectorize(lambda kx, ky, kz : np.squeeze(chi_k([kx, ky, kz]).real))
    interp = interp(kx, ky, kz)
    p.chi_G = np.squeeze(chi_k[Idx(0, 0, 0)].real)    
    p.chi_M = np.squeeze(chi_k[Idx(p.n_k//2, p.n_k//2, 0)].real)    

    if color=='black':
        line = plt.plot(k_plot, interp, '-',c=color)
    else:
        line = plt.plot(k_plot, interp, '-', label=r'$N_{\nu} = $'+'${:d}$'.format(p.nwf),c=color)

    plt.gca().set_xticks(K_plot); plt.gca().set_xticklabels(K_labels)
    plt.grid(True); plt.ylabel(r'$\chi(\mathbf{Q})$')
    plt.legend(loc='best', fontsize=8)
    return line[0].get_color()

import itertools

ws = [0,3,6,] # Selected bosonic frequencies
components = [(aa,bb,cc,dd) for aa,bb,cc,dd in itertools.product([0,1],repeat=4)]

special_K= get_path()[2][1] # M-point on X axis

def comp2int(x):
    return sum(ii*2**nn for nn,ii in enumerate(x[::-1]) )

for w,component in itertools.product(ws,components):
    print(w,component, comp2int( component) )
    if not comp2int(component) in [0,3,6,9,12,15]: continue
    plt.figure(figsize=(3.25*2, 3))

    ps, color = [], None
    for filename in np.sort(glob.glob('data_bse_nwf_*.h5')):
        with HDFArchive(filename, 'r') as a: p = a['p']
        if p.nwf<20: continue
        print(filename)
        subp = [1, 2, 1]
        plt.subplot(*subp); subp[-1] += 1
        color = plot_chi_k(p,p.chi_kw_bse,w,component)
        plt.subplot(*subp); subp[-1] += 1
        plt.plot(1./p.nwf, p.chi_M, 'o', color=color)
        ps.append(p)

    # -- Extrapolation to nwf -> oo
    ps = ParameterCollections(objects=ps)
    x, y = 1./ps.nwf, ps.chi_M
    sidx = np.argsort(x)
    x, y = x[sidx], y[sidx]
    p = np.polyfit(x, y, 1)
    y0 = np.polyval(p, 0)
    X = np.linspace(0, x.max())
    Y = np.polyval(p, X)

    subp = [1, 2, 1]
    plt.subplot(*subp); subp[-1] += 1
    plt.plot(special_K, y0, 'rx')
    plt.subplot(*subp); subp[-1] += 1
    plt.plot(X, Y, '--k', lw=1.0, zorder=-100)
    plt.plot(0, y0, 'rx', label=r'BSE')
    plt.grid(True); plt.legend(loc='best', fontsize=8)
    plt.xlabel(r'$1/N_\nu$'); plt.ylabel(r'$\chi(\mathbf{M})$')
    plt.title(
        r'$\lim_{n_\nu \rightarrow \infty} \, \chi(\mathbf{M}) \approx $' + \
        '${:3.4f}$'.format(y0))
    plt.tight_layout()

    #Now do the same thing for the DBSE
    ps, color = [], None
    for filename in np.sort(glob.glob('data_bse_nwf_*.h5')):
        with HDFArchive(filename, 'r') as a: p = a['p']
        if p.nwf<20: continue
        print(filename)
        subp = [1, 2, 1]
        plt.subplot(*subp); subp[-1] += 1
        color = plot_chi_k(p,p.chi_kw_dbse,w,component,color='black')
        plt.title(f'{component}, w={w}')
        plt.subplot(*subp); subp[-1] += 1
        plt.plot(1./p.nwf, p.chi_M, 'o', color=color)
        ps.append(p)

    # -- Extrapolation to nwf -> oo
    ps = ParameterCollections(objects=ps)
    x, y = 1./ps.nwf, ps.chi_M
    sidx = np.argsort(x)
    x, y = x[sidx], y[sidx]
    p = np.polyfit(x, y, 1)
    y0 = np.polyval(p, 0)
    X = np.linspace(0, x.max())
    Y = np.polyval(p, X)

    subp = [1, 2, 1]
    plt.subplot(*subp); subp[-1] += 1
    plt.plot(special_K, y0, 'rx')
    plt.subplot(*subp); subp[-1] += 1
    plt.plot(X, Y, '--k', lw=1.0, zorder=-100)
    plt.plot(0, y0, 'rx', label=r'DBSE')

    plt.savefig(f'figure_bse_w{w}_comp{comp2int(component)}.svg')
    plt.close()
