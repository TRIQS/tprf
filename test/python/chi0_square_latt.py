
import itertools
import numpy as np

from pytriqs.gf import Gf, MeshImFreq, MeshProduct, MeshBrillouinZone
from pytriqs.lattice import BrillouinZone, BravaisLattice

beta = 10.0
n_k = 50
nw_g = 512
mu = 0.0
nw = 3
nnu = 200

bz = BrillouinZone(BravaisLattice([[1,0],[0,1]]))
bzmesh = MeshBrillouinZone(bz, n_k=n_k)
ek = Gf(mesh=bzmesh, target_shape=[1, 1])

q_list = [q for q in bzmesh]

for idx, k in enumerate(bzmesh):
    #ek[k] = -2*t*(np.cos(k[0]) + np.cos(k[1])) # does not work...
    ek.data[idx] = -2*(np.cos(k[0]) + np.cos(k[1]))

mesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw_g)
bmesh = MeshImFreq(beta=beta, S='Boson', n_max=nw)

iw_list = np.array([ iw for iw in bmesh ])
iw_zero_idx = np.where(iw_list == 0)[0][0]

from pytriqs.applications.tprf.lattice import g0k_from_ek
from pytriqs.applications.tprf.lattice import gr_from_gk
from pytriqs.applications.tprf.lattice import gk_from_gr

g0k = g0k_from_ek(mu=mu, ek=ek, mesh=mesh)
g0r = gr_from_gk(g0k)

from pytriqs.applications.tprf.lattice import chi0r_from_gr_PH
from pytriqs.applications.tprf.lattice import chi0q_from_chi0r
from pytriqs.applications.tprf.lattice import get_at_q
from pytriqs.applications.tprf.lattice import chi0q_sum_nu

chi0r = chi0r_from_gr_PH(nw=nw, nnu=nnu, gr=g0r)
chi0q = chi0q_from_chi0r(chi0r, bz)

if False:
    chi0w0 = np.zeros((len(q_list)), dtype=np.complex)
    for i1, i2 in itertools.product(range(n_k), repeat=2):
        qidx = [i1, i2, 0]
        qidx_lin = bzmesh.index_to_linear(qidx)
        q = np.array(q_list[qidx_lin])/np.pi

        print qidx, q

        chi0 = get_at_q(chi0q, qidx)
        chiw = np.sum(chi0.data, axis=1)
        chiw = np.squeeze(chiw)
        print chiw

        chi0w0[qidx_lin] = chiw
        
chi0w0 = chi0q_sum_nu(chi0q)



from pytriqs.plot.mpl_interface import oplot, oploti, oplotr, plt

# plot zeroth bosonic matsubara freq susceptibility

qx = np.array([q[0] for q in q_list])/np.pi
qy = np.array([q[1] for q in q_list])/np.pi
data = np.squeeze(chi0w0.data[iw_zero_idx])

plt.pcolor(qx.reshape((n_k, n_k)), qy.reshape((n_k, n_k)), data.reshape((n_k, n_k)).real)
plt.colorbar()
plt.axis('equal')
plt.show()

exit()

for q in q_list:
    q = np.array(q)
    q /= np.pi
    plt.plot(q[0], q[1], 'og')

plt.show(); exit()

from pytriqs.plot.mpl_interface import oplot, oploti, oplotr, plt
opt = dict(vmin=-0.01, vmax=0.01, cmap='PuOr')
plt.figure(figsize=(8, 8))

subp = [4, 2, 1]

    

for q in [ [0,0,0], [n_k/2, n_k/2, 0] ]:

    #q = [i1, i2, 0]
    print 'q =', q

    qidx = bzmesh.index_to_linear(q)

    print 'qidx =', qidx
    print 'q_list[qidx] =', q_list[qidx]
    print 'q_list[qidx]/np.pi =', np.array(q_list[qidx])/np.pi
    
    chi0 = get_at_q(chi0q, q)
    chiw = np.sum(chi0.data, axis=1)
    chiw = np.squeeze(chiw)

    print chiw.shape

    #exit()

    plt.subplot(*subp); subp[-1] += 1
    plt.plot(chiw.real)

    plt.subplot(*subp); subp[-1] += 1
    plt.plot(chiw.imag)

    plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(chi0.data[:,:,0,0,0,0].real, **opt)
    plt.colorbar()
    plt.axis('equal')

    plt.subplot(*subp); subp[-1] += 1
    plt.pcolormesh(chi0.data[:,:,0,0,0,0].imag, **opt)
    plt.colorbar()
    plt.axis('equal')

plt.show()
