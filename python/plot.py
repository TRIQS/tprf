
import itertools
import numpy as np

from pytriqs.gf import Idx
from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt

# ----------------------------------------------------------------------
def plot_g2(G2, cplx=None, idx_labels=None, w=Idx(0), opt={}, title=None):

    data = G2[w, :, :].data
    
    if cplx == 're':
        data = data.real
    elif cplx == 'im':
        data = data.imag

    n = data.shape[-1]
    N = n ** 2
    subp = [N, N, 1]

    colorbar_flag = True
    
    import itertools
    for idxs in itertools.product(xrange(n), repeat=4):

        i1, i2, i3, i4 = idxs
        d = data[:, :, i1, i2, i3, i4]
        
        ax = plt.subplot(*subp); subp[-1] += 1

        if idx_labels is not None:
            labels = [ idx_labels[idx] for idx in idxs ]
            sub_title = r'$c^\dagger_{%s} c_{%s} c^\dagger_{%s} c_{%s}$' % tuple(labels)
        else:
            sub_title = str(idxs)

        plt.title(sub_title, fontsize=8)
        
        #plt.pcolormesh(d, **opt)
        if np.max(np.abs(d)) > 1e-10:
            plt.imshow(d, **opt)
            if colorbar_flag:
                if title is not None:
                    plt.title(title)
                plt.colorbar()
                colorbar_flag = False

        ax.set_xticks([])
        ax.set_yticks([])
        plt.axis('equal')


    plt.tight_layout()

# ----------------------------------------------------------------------
def get_g2_opt(G2, cplx=None, cut=1.0):

    data = G2.data
    if cplx == 're':
        data = data.real
    elif cplx == 'im':
        data = data.imag
        
    vmin = np.min(data) * cut
    vmax = np.max(data) * cut
    vm = np.max([np.abs(vmin), np.abs(vmax)])

    opt = dict(vmin=-vm, vmax=vm, cmap=plt.get_cmap('RdBu'))
    
    return opt
