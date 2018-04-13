
import numpy as np
import matplotlib.pyplot as plt
    
# ----------------------------------------------------------------------

from scipy.interpolate import griddata

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive
from pytriqs.gf import MeshBrillouinZone, Idx

# ----------------------------------------------------------------------
def plot_chi(chi):

    bzmesh = chi.mesh.components[1]

    k_vec = np.array([k.value / (2. * np.pi) for k in bzmesh])
    k_vec = k_vec[:, :2]

    # -- Contract Sz Sz response
    Sz = np.diag([+0.5, -0.5])
    chi_SzSz = chi[0, 0, 0, 0].copy()
    chi_SzSz.data[:] = np.einsum('wqabcd,ab,cd->wq', chi.data, Sz, Sz)[:, :, None, None, None, None]

    values = np.squeeze(chi_SzSz[Idx(0), :].data.real)
    print values.shape

    # -- Extend with points beyond the first bz
    k_vec_ext = []
    values_ext = []
    for k_shift in [(0,0), (-1,0), (0, -1), (-1, -1)]:
        k_vec_ext.append( k_vec + np.array(k_shift)[None, :] )
        values_ext.append(values)

    k_vec = np.vstack(k_vec_ext)
    values = np.hstack(values_ext)

    k = np.linspace(-0.75, 0.75, num=400)
    Kx, Ky = np.meshgrid(k, k)
    
    #interp = griddata(k_vec, values, (Kx, Ky), method='nearest')
    #interp = griddata(k_vec, values, (Kx, Ky), method='linear')
    interp = griddata(k_vec, values, (Kx, Ky), method='cubic')

    plt.imshow(
        interp, cmap=plt.get_cmap('magma'),
        extent=(k.min(), k.max(), k.min(), k.max()),
        origin='lower',
        vmin=0, vmax=interp.max(),
        )

    plt.title(r'Lindhardt spin-response $\chi^{00}_{S_z S_z}(\mathbf{q}, \omega=0)$')
    plt.xlabel(r'$\frac{a}{2\pi} \cdot k_x$')
    plt.ylabel(r'$\frac{a}{2\pi} \cdot k_y$')
    plt.colorbar()

# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_ek_and_chi00wq.h5'
    
    with HDFArchive(filename, 'r') as arch:
        chi = arch['chi00wq']
        chi_imtime = arch['chi00wq_imtime']

    plt.figure(figsize=(6, 3))

    subp = [1, 2, 1]

    plt.subplot(*subp); subp[-1] += 1
    plot_chi(chi)

    plt.subplot(*subp); subp[-1] += 1
    plot_chi(chi_imtime)
    
    plt.tight_layout()
    plt.savefig('figure_lindhardt_SzSz_square_latt.pdf')
    plt.show()

