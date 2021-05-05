
""" Non-interacting lattice susceptibility for the square lattice """

import itertools
import numpy as np

from triqs.gf import Gf, MeshImFreq, MeshProduct, MeshBrZone
from triqs.lattice import BrillouinZone, BravaisLattice

from triqs.applications.tprf.lattice import g0k_from_ek
from triqs.applications.tprf.lattice import gr_from_gk

from triqs.applications.tprf.lattice import chi0r_from_gr_PH
from triqs.applications.tprf.lattice import chi0q_from_chi0r
from triqs.applications.tprf.lattice import chi0q_sum_nu

beta = 20.0
n_k = 20
nw_g = 512
mu = 0.0
nw = 10
nnu = 100

# -- BZ-sampling

bz = BrillouinZone(BravaisLattice([[1,0],[0,1]]))
bzmesh = MeshBrZone(bz, n_k=n_k)
q_list = [q for q in bzmesh]

# -- Dispersion

ek = Gf(mesh=bzmesh, target_shape=[1, 1])
for idx, k in enumerate(bzmesh):
    #ek[k] = -2*t*(np.cos(k[0]) + np.cos(k[1])) # does not work...
    ek.data[idx] = -2*(np.cos(k[0]) + np.cos(k[1]))


mesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw_g)
bmesh = MeshImFreq(beta=beta, S='Boson', n_max=nw)

iw_list = np.array([ iw for iw in bmesh ])
iw_zero_idx = np.where(iw_list == 0)[0][0]

# -- Lattice single-particle Green's function

g0k = g0k_from_ek(mu=mu, ek=ek, mesh=mesh)
g0r = gr_from_gk(g0k)

# -- Non-interacting generalized lattice susceptibility

chi0r = chi0r_from_gr_PH(nw=nw, nnu=nnu, gr=g0r)
chi0q = chi0q_from_chi0r(chi0r, bz)

# -- Sum over fermionic Matsubara frequency to get non-int susceptibility

chi0w0 = chi0q_sum_nu(chi0q)

if False: # python implementation
    chi0w0 = np.zeros((len(q_list)), dtype=np.complex)
    for i1, i2 in itertools.product(list(range(n_k)), repeat=2):
        qidx = [i1, i2, 0]
        qidx_lin = bzmesh.index_to_linear(qidx)
        q = np.array(q_list[qidx_lin])/np.pi

        print(qidx, q)

        chi0 = get_at_q(chi0q, qidx)
        chiw = np.sum(chi0.data, axis=1)
        chiw = np.squeeze(chiw)
        print(chiw)

        chi0w0[qidx_lin] = chiw
        
# -- Plot static w=0 susceptibility in k-space

from triqs.plot.mpl_interface import oplot, oploti, oplotr, plt

# plot zeroth bosonic matsubara freq susceptibility

qx = np.array([q[0] for q in q_list])/np.pi
qy = np.array([q[1] for q in q_list])/np.pi
data = np.squeeze(chi0w0.data[iw_zero_idx])

data = -1.0 * data
vmax = np.max(data.real)

plt.title('Square-lattice ' + r'$\chi_0(q, \omega=0)$, $\beta=%2.2f$' % beta)
plt.pcolor(qx.reshape((n_k, n_k)), qy.reshape((n_k, n_k)), data.reshape((n_k, n_k)).real,
           cmap='magma', vmin=0, vmax=vmax)

plt.colorbar()
plt.axis('equal')
plt.xlabel('$q_x/\pi$')
plt.ylabel('$q_y/\pi$')
plt.tight_layout()
plt.savefig('figure_chi0q_w0_square_latt.pdf')

# -- Plot Gamma and M point dynamic susceptibilities

opt = dict(vmin=-0.01, vmax=0.01, cmap='PuOr')
plt.figure(figsize=(12, 8))

subp = [3, 4, 1]    

for q in [ [0,0,0], [n_k/2, n_k/2, 0], [n_k/2, 0, 0] ]:

    qidx = bzmesh.index_to_linear(q)

    print('-'*72)
    print('q =', q)
    print('qidx =', qidx)
    print('q_list[qidx] =', q_list[qidx])
    print('q_list[qidx]/np.pi =', np.array(q_list[qidx])/np.pi)
    
    data = np.squeeze(chi0w0.data[:, qidx])
    print(data.shape)
    
    plt.subplot(*subp); subp[-1] += 1
    plt.title(r'$q = \pi \times$ %s' % str(np.array(q_list[qidx])/np.pi))
    plt.plot(data.real)
    plt.ylabel(r'Re[$\chi_0(i\omega)$]')
    plt.ylim([-vmax, 0.1*vmax])

    plt.subplot(*subp); subp[-1] += 1
    plt.plot(data.imag)
    plt.ylabel(r'Im[$\chi_0(i\omega)$]')

    plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(chi0q.data[:,:,qidx,0,0,0,0].real, **opt)
    plt.colorbar()
    plt.axis('equal')
    plt.title(r'Re[$\chi_0(i\omega, i\nu)$]')
    plt.xlabel(r'$\nu$')
    plt.ylabel(r'$\omega$')

    plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(chi0q.data[:,:,qidx,0,0,0,0].imag, **opt)
    plt.colorbar()
    plt.axis('equal')
    plt.title(r'Im[$\chi_0(i\omega, i\nu)$]')
    plt.xlabel(r'$\nu$')
    plt.ylabel(r'$\omega$')

plt.tight_layout()
plt.savefig('figure_chi0q_w_square_latt.pdf')

plt.show()
