import itertools
import types

import numpy as np

import matplotlib as mpl
from pytriqs.plot.mpl_interface import plt

from triqs_tprf.lattice_utils import k_space_path

hs_to_k = {
                'G' : np.array([0.0, 0.0, 0.0]), 
                'X' : np.array([0.5, 0.0, 0.0]),
                'Y' : np.array([0.0, 0.5, 0.0]), 
                'Z' : np.array([0.0, 0.0, 0.5]),
                'M' : np.array([0.5, 0.5, 0.0]), 
                'R' : np.array([0.5, 0.5, 0.5]), 
                        }

hs_to_latex = {
                'G' : r'$\Gamma$',
                'X' : r'$X$',
                'Y' : r'$Y$',
                'Z' : r'$Z$',
                'M' : r'$M$',
                'R' : r'$R$',
                }

def bsplot(obj, path, *opt_list, **opt_dict):
    """
    Plot stuff like bs lol
    """
    __bsplot_impl(plt, obj, path, plt.xticks, plt.xticks, *opt_list, **opt_dict)
        

def __bsplot_impl(top, obj, path, xticks_fct, xticklabels_fct, *opt_list, **opt_dict):

    hs_points = path.split('-')
    hs_labels = [hs_to_latex[hs_point] for hs_point in hs_points]
    hs_k = [hs_to_k[hs_point] for hs_point in hs_points]

    k_paths = zip(hs_k, hs_k[1:])
    k_vecs, k_plot, K_plot = k_space_path(k_paths, bz=obj.mesh.domain)
    kx, ky, kz = k_vecs.T
    
    plt_fct = getattr(top, 'plot')

    get_gf_on_path = np.vectorize(lambda kx, ky, kz : obj([kx, ky, kz]).real)    
    gf_on_path = get_gf_on_path(kx, ky, kz)

    plt_fct(k_plot, gf_on_path, *opt_list, **opt_dict)
    
    if isinstance(top, types.ModuleType):
        xticks_fct(K_plot, hs_labels)
    else:
        xticks_fct(K_plot)
        xticklabels_fct(hs_labels)

mpl.axes.Axes.bsplot = lambda self, obj, path, *opt_list, **opt_dict: __bsplot_impl(self, obj, path, self.set_xticks, self.set_xticklabels, *opt_list, **opt_dict)




