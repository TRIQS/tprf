from common import *
from triqs.plot.mpl_interface import plt,oplot
import triqs.gf
from h5 import HDFArchive
import itertools
# ----------------------------------------------------------------------

from h5 import HDFArchive
from triqs.gf import MeshBrillouinZone, Idx, Gf, MeshProduct
from triqs.lattice import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

from triqs.lattice.utils import k_space_path

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections

# ----------------------------------------------------------------------

def load_h5(filename):
    print(f'--> Loading: {filename}')
    with HDFArchive(filename, 'r') as a:
        p = a['p']
    return p

p=load_h5('data_vertices.h5')

font = {'size'   : 8}

import matplotlib
matplotlib.rc('font', **font)


norbs = 2
beta=5 

def plot_vertex(norbs=norbs,W=0):
    print('vertex')
    scale= p.U/(beta*beta)*1 # Note normalization with beta^2

    fig,axs=plt.subplots(nrows= norbs**2, ncols=norbs**2, figsize=(10,10))
    for aa,bb,cc,dd in itertools.product( range(norbs),repeat=4):
        index = (aa*norbs+bb, cc*norbs+dd)
        z_plot=axs[index].imshow(p.Gamma_wnn[Idx(W),:,:].data[:,:,aa,bb,cc,dd].real,
            vmax=scale,vmin=-scale,cmap=plt.get_cmap('terrain'),
            )
        plt.colorbar(z_plot,ax=axs[index])

    plt.savefig('fig_Gamma_array.svg')

    fig,axs=plt.subplots(nrows= norbs**2, ncols=norbs**2, figsize=(10,10))
    for aa,bb,cc,dd in itertools.product( range(norbs),repeat=4):
        index = (aa*norbs+bb, cc*norbs+dd)
        z_plot=axs[index].imshow(p.F_wnn[Idx(W),:,:].data[:,:,aa,bb,cc,dd].real,
            vmax=scale,vmin=-scale,cmap=plt.get_cmap('terrain'),
            )
        plt.colorbar(z_plot,ax=axs[index])

    plt.savefig('fig_F_array.svg')

    return

plot_vertex(W=1)
