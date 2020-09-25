# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019, S. Käser
# Author: S. Käser
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

import itertools
import types

import numpy as np

import matplotlib as mpl
from triqs.plot.mpl_interface import plt

from triqs_tprf.lattice_utils import k_space_path

from scipy.stats import gaussian_kde

# ========== Bandstructure ========== 

def bsplot(obj, path, *opt_list, **opt_dict):
    """
    Plot Gf objects like bandstructure
    """
    __bsplot_impl(plt, obj, path, *opt_list, **opt_dict)
        
def __bsplot_impl(top, obj, path, *opt_list, **opt_dict):

    hs_labels, hs_k = list(zip(*path))

    k_paths = list(zip(hs_k, hs_k[1:]))
    k_vecs, k_plot, K_plot = k_space_path(k_paths, bz=obj.mesh.domain)
    kx, ky, kz = k_vecs.T

    # -- If top is plt do now access the currently used Axes object
    if isinstance(top, types.ModuleType):
        top = top.gca()

    get_gf_on_path = np.vectorize(lambda kx, ky, kz : obj([kx, ky, kz]).real)    
    gf_on_path = get_gf_on_path(kx, ky, kz)

    top.plot(k_plot, gf_on_path, *opt_list, **opt_dict)

    top.set_xticks(K_plot)
    top.set_xticklabels(hs_labels)

    # -- Make x-spine cut off at starting and ending high symmetry point
    top.spines['bottom'].set_bounds(top.get_xticks()[0], top.get_xticks()[-1])

    # -- Make y-spine cut off at highest and lowest value 
    lower_limit = np.min([line.get_ydata() for line in top.get_lines()])
    upper_limit = np.max([line.get_ydata() for line in top.get_lines()])

    top.spines['left'].set_bounds(lower_limit, upper_limit)
    top.set_yticks([lower_limit, upper_limit])

mpl.axes.Axes.bsplot = lambda self, obj, path, *opt_list, **opt_dict : \
        __bsplot_impl(self, obj, path, *opt_list, **opt_dict)

# ========== DOS ========== 

def dosplot(obj, *opt_list, **opt_dict):
    """Plot density of states for dispersion relation objects
    """

    __dosplot_impl(plt, obj, *opt_list, **opt_dict)

def __dosplot_impl(top, obj, *opt_list, **opt_dict):

    lower_limit = np.min(obj.data.real)
    upper_limit = np.max(obj.data.real)

    dos = gaussian_kde(obj.data.real)
    xs = np.linspace(lower_limit, upper_limit, 500)

    dos.covariance_factor = lambda : .1
    dos._compute_covariance()

    # -- If top is plt do now access the currently used Axes object
    if isinstance(top, types.ModuleType):
        top = top.gca()

    top.plot(dos(xs).real, xs, *opt_list, **opt_dict)
    top.fill_betweenx(xs, dos(xs).real, [0]*len(xs), alpha=0.25, *opt_list, **opt_dict)

    # -- Make y-spine cut off at highest and lowest value 
    lower_limit = np.min([line.get_ydata() for line in top.get_lines()])
    upper_limit = np.max([line.get_ydata() for line in top.get_lines()])

    top.spines['left'].set_bounds(lower_limit, upper_limit)
    top.set_yticks([lower_limit, upper_limit])

    # -- No x-ticks
    top.set_xticks([])

mpl.axes.Axes.dosplot = lambda self, obj, *opt_list, **opt_dict : \
        __dosplot_impl(self, obj, *opt_list, **opt_dict)
