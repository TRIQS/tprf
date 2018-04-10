
import numpy as np
import matplotlib.pyplot as plt
    
# ----------------------------------------------------------------------

from scipy.interpolate import griddata
from skimage.measure import find_contours

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive
from pytriqs.gf import MeshBrillouinZone

# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_ek_and_chi00wq.h5'
    
    with HDFArchive(filename, 'r') as arch:
        ek = arch['ek']

    bzmesh = ek.mesh

    k_vec = np.array([k.value / (2. * np.pi) for k in bzmesh])
    k_vec = k_vec[:, :2]
    
    values = ek.data[:, 0, 0].real

    # -- Extend with points beyond the first bz
    k_vec_ext = []
    values_ext = []
    for k_shift in [(0,0), (-1,0), (0, -1), (-1, -1)]:
        k_vec_ext.append( k_vec + np.array(k_shift)[None, :] )
        values_ext.append(values)

    k_vec = np.vstack(k_vec_ext)
    values = np.hstack(values_ext)

    #k = np.linspace(k_vec.min(), k_vec.max(), num=400)
    k = np.linspace(-0.5, 0.5, num=400)
    Kx, Ky = np.meshgrid(k, k)
    
    ek_interp = griddata(k_vec, values, (Kx, Ky), method='cubic')

    fermi_surfs = find_contours(ek_interp, 0.0)

    plt.imshow(
        ek_interp, cmap=plt.get_cmap('RdBu'),
        extent=(k.min(), k.max(), k.min(), k.max()),
        origin='lower',
        )

    for fidx, fs in enumerate(fermi_surfs):
        fs /= ek_interp.shape[0]
        fs += np.array([k.min()]*2)[None, :]
        plt.plot(fs[:, 0], fs[:, 1], '-k', lw=0.5)

    plt.xlabel(r'$\frac{a}{2\pi} \cdot k_x$')
    plt.ylabel(r'$\frac{a}{2\pi} \cdot k_y$')

    plt.tight_layout()
    plt.savefig('figure_fermi_surf_square_latt.pdf')
    plt.show()
