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

from common import *
    
filename = './data/data_sc.h5'
print(f'--> Loading: {filename}')

with HDFArchive(filename, 'r') as a:
    ps = ParameterCollections(a['ps'])
p = ps.objects[-1]

filename_ref = './data/data_g.h5'
print(f'--> Loading: {filename_ref}')
with HDFArchive(filename_ref, 'r') as a:
    p_ref_g = a['p']

# -----------------------------------

from triqs.plot.mpl_interface import oplot, oplotr, oploti, plt

plt.figure(figsize=(3.25*4, 10))
subp = [5, 2, 1]

if True:
    plt.subplot(*subp); subp[-1] += 1

    #G = np.array([0.0, 0.0, 0.0])
    #M = np.array([0.5, 0.0, 0.0])
    #K = np.array([1.0/3, 1.0/3, 0.0])

    #path = [(G,M), (M,K), (K,G)]
    #labels = [r'$\Gamma$', '$M$', '$K$', r'$\Gamma$']

    Gamma, X, M = [0.00, 0.00, 0.00], [ 0.00, 0.00, 0.50], [-0.25, 0.25, 0.25]
    path = [(Gamma, X), (X, M), (M, Gamma),]
    labels = [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$', ]
    
    from triqs.lattice.utils import k_space_path

    k_vecs, k_plot, k_ticks = k_space_path(
        path, num=32, bz=p.e_k.mesh.bz, relative_coordinates=False)

    e_k_interp = np.vectorize(lambda k : np.linalg.eigvalsh(p.e_k(k)).real, signature='(n)->(m)')

    plt.plot(k_plot, e_k_interp(k_vecs))
    plt.xticks(k_ticks, labels=labels)
    plt.ylabel(r'$\epsilon(\mathbf{k})$')
    plt.grid(True)

plt.subplot(*subp); subp[-1] += 1
plt.plot(ps.iter, ps.dG, 's-')
plt.plot(ps.iter, ps.iter * 0 + p.G_tol, '-r', lw=1)
plt.ylabel('$\max | \Delta G_l |$')
plt.xlabel('Iteration')
plt.semilogy([], [])

plt.subplot(*subp); subp[-1] += 1
plt.plot(ps.iter, ps.mu, 's-')
plt.ylabel('$\mu$')
plt.xlabel('Iteration')

plt.subplot(*subp); subp[-1] += 1
plt.plot(ps.iter, ps.N, 's-')
if hasattr(p, 'N_target'):
    plt.plot(ps.iter, ps.iter * 0 + p.N_target, '-k', lw=1)
    plt.plot(ps.iter, ps.iter * 0 + p.N_target + p.N_tol, '-r', lw=1)
    plt.plot(ps.iter, ps.iter * 0 + p.N_target - p.N_tol, '-r', lw=1)
plt.ylabel('$N$')
plt.xlabel('Iteration')

plt.subplot(*subp); subp[-1] += 1
rho = np.diagonal(ps.rho.real, axis1=1, axis2=2)
plt.plot(ps.iter, rho[:, 0], 'o-', label='n_xy')
plt.plot(ps.iter, rho[:, 2], 'o-', label='n_xz, n_yz')
plt.legend(loc='best')
plt.xlabel('Iteration')

plt.subplot(*subp); subp[-1] += 1
g_xy_dlr = p.G_tau['up_2']
g_xz_dlr = p.G_tau['up_1']
g = p.G_tau_raw
g_xy = (g['up_2'] + g['do_2']) / 2
g_xz = (g['up_1'] + g['do_1'] + g['up_0'] + g['do_0']) / 4
oplotr(-g_xy_dlr, label=f'DLR xy')
oplotr(-g_xz_dlr, label=f'DLR xz,yz')
oplotr(-g_xy, '-k', alpha=0.25, label=f'raw xy')
oplotr(-g_xz, '-k', alpha=0.25, label=f'raw xz,yz')
plt.legend(loc='best', fontsize=8)
plt.ylabel(r'$G(\tau)$')

plt.subplot(*subp); subp[-1] += 1
oplotr(-g_xy_dlr + g_xy, label=f'DLR dG xy')
oplotr(-g_xz_dlr + g_xz, label=f'DLR dG xz,yz')

plt.subplot(*subp); subp[-1] += 1
oplotr(p.sigma_w[0,0], '.-', label='S_xy', lw=0.5)
oplotr(p.sigma_w[1,1], '.-', label='S_xz', lw=0.5)
oplotr(p.sigma_w[2,2], '.-', label='S_yz', lw=0.5)
plt.ylabel(r'$\Sigma(i\omega_n)$')
plt.xlim([0, 10])

plt.subplot(*subp); subp[-1] += 1
oploti(p.sigma_w[0,0], '.-', label='S_xy', lw=0.5)
oploti(p.sigma_w[1,1], '.-', label='S_xz', lw=0.5)
oploti(p.sigma_w[2,2], '.-', label='S_yz', lw=0.5)
plt.ylabel(r'$\Sigma(i\omega_n)$')
plt.xlim([0, 10])
plt.ylim(top=0)

plt.tight_layout()
plt.savefig('figure_sc.pdf')
plt.show()
