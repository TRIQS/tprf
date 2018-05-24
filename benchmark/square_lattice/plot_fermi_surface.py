
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

    filename = 'data_e_k_and_chi00_wk.h5'
    
    with HDFArchive(filename, 'r') as arch:
        e_k = arch['e_k']
            
    k = np.linspace(-0.5, 0.5, num=200) * 2. * np.pi
    Kx, Ky = np.meshgrid(k, k)
    
    e_k_interp = np.vectorize(lambda kx, ky : e_k([kx, ky, 0])[0,0].real)    
    e_k_interp = e_k_interp(Kx, Ky)

    plt.imshow(
        e_k_interp, cmap=plt.get_cmap('RdBu'),
        extent=(k.min(), k.max(), k.min(), k.max()),
        origin='lower',
        )
    plt.colorbar()
    
    plt.contour(Kx, Ky, e_k_interp, levels=[0])
    
    plt.xlabel(r'$\frac{a}{2\pi} \cdot k_x$')
    plt.ylabel(r'$\frac{a}{2\pi} \cdot k_y$')

    plt.tight_layout()
    plt.savefig('figure_fermi_surf_square_latt.pdf')
    plt.show()
