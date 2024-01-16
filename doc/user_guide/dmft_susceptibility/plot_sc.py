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

from triqs.gf import make_gf_imtime

with HDFArchive('./data/data_sc.h5', 'r') as a: ps = a['ps']
p = ps.objects[-1]

#with HDFArchive('./data/data_non_int.h5', 'r') as a: p_non_int = a['p']

from fitdlr import triqs_driver
dd = triqs_driver(p.g_c.mesh)

tmp = p.sigma_w.copy()
tmp.data[:] -= p.S0[None, ...]
p.sigma_c_fit = make_gf_dlr(tmp)
p.sigma_w_dense = make_gf_imfreq(p.sigma_c_fit, n_iw=100) + p.S0

iwn = np.array([complex(x) for x in p.Sigma_w.mesh])

#S_high_freq = p.Sigma_moments['0'][0] + p.Sigma_moments['0'][1]/iwn[:, None, None]
#G_high_freq = p.G_moments['0'][0] + p.G_moments['0'][1]/iwn[:, None, None] + p.G_moments['0'][2]/iwn[:, None, None]**2 

#for key, value in p.G_moments.items():
#    print(f'key = {key}\n{value}')

#exit()

S_high_freq = np.array([ p.Sigma_moments[b][0][0,0] + p.Sigma_moments[b][1][0,0]/iwn for b in p.spin_orb_names ])
G_high_freq = np.array([ p.G_moments[b][0][0,0] + p.G_moments[b][1][0,0]/iwn + p.G_moments[b][2][0,0]/iwn**2 for b in p.spin_orb_names ])

#print(S_high_freq.shape)
#print(G_high_freq.shape)
#exit()
#exit()
#exit()

def lims(y=None):
    if y is not None:
        plt.ylim(y)
    plt.xlim([-250, 250])

#plt.figure(figsize=(3.25*2, 5))
plt.figure(figsize=(10, 10))
subp = [5, 2, 1]

plt.subplot(*subp); subp[-1] += 1
plt.plot(ps.iter, ps.dG, 's-')
plt.ylabel('$\max | \Delta G(\tau) |$')
plt.xlabel('Iteration')
plt.semilogy([], [])

plt.subplot(*subp); subp[-1] += 1
colors = [ plt.plot([], [])[0].get_color() for x in range(2) ]
keys = ['up_0', 'do_0']
oplotr(p.G_tau_raw, alpha=0.25, label=None)
for i, c in enumerate(colors):
    oplotr(p.G_tau[keys[i]], label=None, color=c)
    #oplotr(p_non_int.g_tau[i,i], ':', lw=5, label=None, color=c)
plt.legend(loc='best', fontsize=8)
plt.ylabel(r'$G(\tau)$')

plt.subplot(*subp); subp[-1] += 1
for i, j in itertools.product(range(2), repeat=2):
    plt.plot(dd.w_x, p.g_c[i, j].data.flatten().real, '.b')
    plt.plot(dd.w_x, -p.g_c[i, j].data.flatten().real, '.r')
plt.semilogy([], [])

plt.subplot(*subp); subp[-1] += 1
oplotr(p.G_tau_raw - p.G_tau, alpha=0.25, label=None)
oplotr(-p.G_tau_raw + p.G_tau, alpha=0.25, label=None)
plt.grid(True)
plt.ylabel('DLR fit difference')
plt.xlabel(r'$\tau$')

plt.subplot(*subp); subp[-1] += 1
oplotr(p.G0_w)
lims()

plt.subplot(*subp); subp[-1] += 1
oploti(p.G0_w)
lims()

plt.subplot(*subp); subp[-1] += 1
for i in range(p.num_spin_orb):
    plt.plot(iwn.imag, G_high_freq.real[i], ':')
oplotr(p.G_w, '+', label=None)
for i in range(2):
    oplotr(p.g_w[i, i], '.', label=None)
    oplotr(p.g_w_linear[i, i], label=None)
lims()

plt.subplot(*subp); subp[-1] += 1
for i in range(p.num_spin_orb):
    plt.plot(iwn.imag, G_high_freq.imag[i], ':')
oploti(p.G_w, '+', label=None)
for i in range(2):
    oploti(p.g_w[i, i], '.', label=None)
    oploti(p.g_w_linear[i, i], label=None)
lims()

plt.subplot(*subp); subp[-1] += 1
for i in range(p.num_spin_orb):
    plt.plot(iwn.imag, S_high_freq.real[i], ':')
    oplotr(p.sigma_w[i, i], label=None)
    oplotr(p.sigma_w_dense[i, i], label=None)

oplotr(p.Sigma_w, '-', alpha=0.5, label=None)
plt.ylabel(r'$\Sigma(i\omega_n)$')
lims(y=[-4, 10])

plt.subplot(*subp); subp[-1] += 1
for i in range(p.num_spin_orb):
    plt.plot(iwn.imag, S_high_freq.imag[i], ':')
    oploti(p.sigma_w[i, i], label=None)
    oploti(p.sigma_w_dense[i, i], label=None)

oploti(p.Sigma_w, '-', alpha=0.5, label=None)
plt.ylabel(r'$\Sigma(i\omega_n)$')
lims(y=[-5, 5])

plt.tight_layout()
plt.savefig('figure_sc.svg')
plt.show()
