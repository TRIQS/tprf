
""" Test of higher order fourier transform """

# ----------------------------------------------------------------------

import itertools
import numpy as np

from triqs.gf import Gf, MeshImFreq, MeshImTime, MeshProduct

# ----------------------------------------------------------------------
def first_index(mesh):
    if hasattr(mesh, 'first_index'):
        return mesh.first_index()
    else:
        return 0

def mesh_product_iterator(mesh):
    return itertools.product(
        *[enumerate(mesh, start=first_index(mesh))
          for mesh in mesh.components])

def mesh_product_iterator_numpy(mesh):
    return [list(map(np.array, list(zip(*list(x)))))
            for x in mesh_product_iterator(chi4_tau.mesh)]

# ----------------------------------------------------------------------
if __name__ == '__main__':

    nw = 20
    nt = 45
    beta = 2.0

    imtime = MeshImTime(beta, 'Fermion', nt)
    imfreq = MeshImFreq(beta, 'Fermion', nw)

    imtime3 = MeshProduct(imtime, imtime, imtime)
    imfreq3 = MeshProduct(imfreq, imfreq, imfreq)

    chi4_tau = Gf(name=r'$g(\tau)$', mesh=imtime3, target_shape=[1, 1, 1, 1])

    print(dir(chi4_tau.indices))
    for i in chi4_tau.indices:
        print(i)
    exit()
    
    # -- Smooth anti-periodic function
    w = 0.5
    e1, e2, e3 = w, w, w
    E = np.array([e1, e2, e3])
    for idx, tau in mesh_product_iterator_numpy(chi4_tau.mesh):
        chi4_tau[idx.tolist()][:] = np.sum(np.cos(np.pi*E*tau))

    # -- Test fourier

    from triqs.applications.susceptibility.fourier import chi4_iw_from_tau
    chi4_iw = chi4_iw_from_tau(chi4_tau, nw)

    from triqs.applications.susceptibility.fourier import chi4_tau_from_iw
    chi4_tau_ref = chi4_tau_from_iw(chi4_iw, nt)
    
    np.testing.assert_array_almost_equal(chi4_tau_ref.data, chi4_tau.data)

    # ------------------------------------------------------------------
    if False:
        from triqs.plot.mpl_interface import oplot, plt

        subp = [2, 2, 1]
        plt.subplot(*subp); subp[-1] += 1
        plt.title('chi4_tau')
        cut = chi4_tau[[0, all, all]]
        plt.pcolormesh(cut.data[:,:,0,0,0,0].real)
        plt.axis('square')

        plt.subplot(*subp); subp[-1] += 1
        plt.title('chi4_tau_ref')
        cut = chi4_tau_ref[[0, all, all]]
        plt.pcolormesh(cut.data[:,:,0,0,0,0].real)
        plt.axis('square')
        
        cidx = 5 # cut idx
        cut = chi4_iw[[cidx, all, all]]
         
        plt.subplot(*subp); subp[-1] += 1
        plt.title('chi4_iw Re')
        plt.pcolormesh(cut.data[:,:,0,0,0,0].real)
        plt.colorbar()
        plt.axis('square')

        plt.subplot(*subp); subp[-1] += 1
        plt.title('chi4_iw Im')
        plt.pcolormesh(cut.data[:,:,0,0,0,0].imag)
        plt.colorbar()
        plt.axis('square')

        plt.tight_layout()
        plt.show()
