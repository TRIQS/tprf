
import os
import copy
import glob
import itertools
import numpy as np
import matplotlib.pyplot as plt

#import ase.units as units

# ----------------------------------------------------------------------

from h5 import HDFArchive
from triqs.gf import MeshBrillouinZone, Idx, Gf, MeshProduct
from triqs.lattice import BrillouinZone, BravaisLattice

# ----------------------------------------------------------------------

from triqs.lattice.utils import k_space_path

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections

# ----------------------------------------------------------------------

from common import *

# ----------------------------------------------------------------------
def do_post_processing(p):

    chi_wk = p.chi_wk
    bzmesh = p.chi_wk.mesh.components[1]
        
    #p.n_w = p.nwf
    p.beta = p.chi_wk.mesh[0].beta
    
    # -- Plot band structure alon a high symmetry path in the Brillouin Zone

    Gamma, X, M = [0.00, 0.00, 0.00], [ 0.00, 0.00, 0.50], [-0.25, 0.25, 0.25]
    paths = [(Gamma, X), (X, M), (M, Gamma),]
    labels = [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$', ]

    k_vecs, k_plot, K_plot = k_space_path(paths, bz=bzmesh.bz, num=10000, relative_coordinates=False)
    k_vecs_abs = k_vecs
    
    # ------------------------------------------------------------------
    # -- Spin-spin response

    def get_O1O2(chi_wk, O1, O2, label):

        oo = ParameterCollection()        
        oo.label = label
        
        chi_O1O2 = chi_wk[Idx(0), :][0, 0, 0, 0].copy()        
        chi_O1O2.data[:] = np.einsum(
            'ab,qabcd,cd->q', O1, chi_wk[Idx(0), :].data, O2)

        chi_interp = np.zeros(
            [k_vecs.shape[0]] + list(chi_O1O2.target_shape), dtype=complex)

        for kidx, (kx, ky, kz) in enumerate(k_vecs_abs):
            chi_interp[kidx] = chi_O1O2((kx, ky, kz))

        s_chi_O1O2 = lambda k : chi_O1O2( np.array(k) @ bzmesh.bz.units )
                    
        oo.chi_G = s_chi_O1O2(Gamma)
        oo.chi_M = s_chi_O1O2(M)
        oo.chi_X = s_chi_O1O2(X)
        
        oo.chi_interp = chi_interp.real

        return oo
    
    # ------------------------------------------------------------------

    kx, ky, kz = k_vecs_abs.T
    
    o = ParameterCollection()
    
    o.n_k = p.n_k
    #o.n_w = p.n_w
    o.nwf = p.nwf
    o.label = p.label

    o.beta = p.beta
    kB = 8.61733033722e-05
    o.T = 1./p.beta/kB

    sz = 0.5*np.diag([1., -1.])
    Sz = np.kron(sz, np.eye(3))

    o.SzSz = get_O1O2(p.chi_wk, Sz, Sz, r'$\langle S_z S_z \rangle$')

    o.k_vecs = k_vecs
    o.k_plot = k_plot
    o.K_plot = K_plot

    o.label = p.label
    o.labels = labels

    return o

# ----------------------------------------------------------------------
def read_postprocess_and_write_depr(filename):

    d = ParameterCollections()
    
    print('--> Loading:', filename)

    def swap_kw_to_wk(g_kw):
        km, wm = g_kw.mesh.components
        g_wk = Gf(mesh=MeshProduct(wm, km), target_shape=g_kw.target_shape)
        g_wk.data[:] = g_kw.data.swapaxes(0, 1)
        return g_wk
    
    with HDFArchive(filename, 'r') as arch:

        p = copy.deepcopy(arch['p'])
        p.label = 'DFT+DMFT+BSE'
        p.chi_wk = swap_kw_to_wk(p.chi_kw_bse)

        # --
        
        pd = copy.deepcopy(arch['p'])
        pd.label = 'DFT+DMFT+D-BSE'
        pd.chi_wk = swap_kw_to_wk(p.chi_kw_dbse)

        # --
        
        p0 = copy.deepcopy(arch['p'])
        p0.label = 'DFT+DMFT'
        p0.chi_wk = swap_kw_to_wk(p.chi0_kw)

    o = do_post_processing(p)        

    out_filename = filename.replace('.h5', '_dmft_bse_interp.h5')
    print('--> Storing:', out_filename)
    with HDFArchive(out_filename, 'w') as arch:
        arch['d'] = o

    od = do_post_processing(pd)        

    out_filename = filename.replace('.h5', '_dmft_dbse_interp.h5')
    print('--> Storing:', out_filename)
    with HDFArchive(out_filename, 'w') as arch:
        arch['d'] = od
        
    o0 = do_post_processing(p0) 

    out_filename = filename.replace('.h5', '_dmft_interp.h5')
    print('--> Storing:', out_filename)
    with HDFArchive(out_filename, 'w') as arch:
        arch['d'] = o0

# ----------------------------------------------------------------------
def read_postprocess_and_write(filename):

    print('--> Assessing:', filename)

    out_filename = filename.replace('.h5', '_interp.h5')

    if os.path.isfile(out_filename):
        print('--> Interpolated data exist, skipping..')
        #return

    print('--> Loading:', filename)

    with HDFArchive(filename, 'r') as arch:
        p = arch['p']
        
    o = do_post_processing(p)        

    print('--> Storing:', out_filename)
    with HDFArchive(out_filename, 'w') as arch:
        arch['d'] = o
        
# ----------------------------------------------------------------------
def plot(filename, style='-', color=None, static_response=False, label=False, skip_list=[], marker='o'):
    
    print('--> Loading:', filename)
    with HDFArchive(filename, 'r') as arch: p = arch['d']

    #if p.n_w in skip_list:
    #    return color, p.n_w, p

    if hasattr(p, 'n_w'):
        p.nwf = p.n_w
    
    if p.nwf in skip_list:
        return color, p.nwf, p

    ax = plt.gca()
    
    # -- Approx static susceptiblities
    chi_G = 5.393
    chi_M = 4.086
    #chi_33 = 9.489
    #chi_33 = 8.255
    chi_33 = 8.099620682819879 * 3./2.

    #K_G = p.K_plot[0]
    #K_M = p.K_plot[1]
    #K_X = p.K_plot[2]
    #K_G2 = p.K_plot[3]    
    #K_33 = K_X + 1./3. * (K_G2 - K_X)

    if len(p.K_plot) < 9:
        K_G = p.K_plot[0]
        K_X = p.K_plot[1]
        K_M = p.K_plot[2]
        K_33 = K_G + 2./3. * (K_X - K_G)
    else:
        K_G = p.K_plot[1]
        K_X = p.K_plot[2]
        K_M = p.K_plot[3]
        #K_G2 = p.K_plot[4]    
        K_33 = K_G + 2./3. * (K_X - K_G)

    if static_response:
        plt.plot([K_G, K_X, K_33], [chi_G, chi_M, chi_33], 'ro',
                 alpha=0.5, zorder=100,
                 #label=r'$\chi_{S_z S_z}$ SC-DMFT'
                 label=r'$\chi$ SC-DMFT'
                 )
    
    comp = [p.SzSz]
    for oo in comp:

        #print(oo.chi_loc)
        #print(oo.chi_interp.shape)

        if color is not None:
            opts = dict(color=color)
        else:
            opts = dict()
        
        if label:
            #label = r'%s $\chi_{loc} = %2.2f$' % (oo.label, oo.chi_loc.real)
            #label = r'%s, $N_\nu = %d$, %s' % (oo.label, p.n_w, p.label)
            #label = r'$N_\nu = %d$' % (p.n_w)
            label = r'$N_\nu = %d$' % (p.nwf)
        else:
            label = None
        lines = plt.plot(p.k_plot, oo.chi_interp, style, alpha=0.75, label=label, **opts)
        color = lines[0].get_color()

        for K, chi in [(K_X, oo.chi_X)]:
            plt.plot(K, chi, marker, color=color, alpha=0.75, markersize=4)

        #for K, chi in [(K_G, oo.chi_G), (K_M, oo.chi_M), (K_X, oo.chi_X), (K_33, oo.chi_IC)]:
        #    plt.plot(K, chi, '.', color=color)
                    
    plt.grid()
    ax.set_xticks(p.K_plot)
    plt.xlim([p.K_plot.min(), p.K_plot.max()])
    #ax.set_xticklabels([r'$\Gamma$',r'$M$',r'$X$',r'$\Gamma$', r'$Z$'])
    ax.set_xticklabels(p.labels)

    if len(p.K_plot) > 9:
        plt.text(p.K_plot[6], -0.1, r'$Y_0$', ha='left', va='bottom')
        plt.text(p.K_plot[9], -0.1, r'$M$', ha='left', va='bottom')
    
    plt.legend(loc='upper right', ncol=3)
    
    #plt.ylim([-0.1, 2])
    #plt.ylim([-0.1, p.SzSz.chi_interp.max() * 1.4])
    #plt.ylim(bottom=0)
    #plt.savefig('figure_chi00_bandpath_%s.pdf' % p.label)
    #plt.ylabel(r'$\chi_{S_z S_z}(\mathbf{q}, \omega = 0)$')
    plt.ylabel(r'$\chi_{S_z S_z}(\mathbf{q})$')

    return color, p.nwf, p
        

if __name__ == '__main__':

    nk = '048'

    nwfs = [f'{i:03d}' for i in np.arange(5, 45, 5)]
    skip_list = []
        
    if False:
        for nwf in ['004']:
            filename = f'./data/data_bse_nwf_{nwf}_nk_{nk}.h5'
            read_postprocess_and_write_depr(filename)

    plt.figure(figsize=(3.25*2, 3.75*1.5))

    from matplotlib.gridspec import GridSpec
    gs = GridSpec(
        2, 1,
        width_ratios=[1],
        height_ratios=[1, 0.5],
        wspace=0.0, hspace=0.2,
        bottom=0.10, top=0.98,
        left=0.08, right=0.98,
        )
    
    plt.subplot(gs[0, 0])

    plt.plot([], [], 'o-', color='gray', label=r'DMFT+DBSE', alpha=0.75, markersize=4)
    plt.plot([], [], 's--', color='gray', label=r'DMFT+BSE', alpha=0.75, markersize=4)

    kind = 'dmft_dbse'
    style = '-'

    colors = dict()
    psd = []
    for nwf in nwfs:
        c, n, p = plot(
            f'./data/data_bse_nwf_{nwf}_nk_{nk}_{kind}_interp.h5',
            style=style, static_response=False, label=True, skip_list=skip_list)
        colors[n] = c
        psd.append(p)

    kind = 'dmft_bse'
    style = '--'
    ps = []
    for nwf in nwfs:
        n = int(nwf.split('_')[0])
        _, _, p = plot(
            f'./data/data_bse_nwf_{nwf}_nk_{nk}_{kind}_interp.h5',
            style=style, static_response=False, color=colors[n], skip_list=skip_list, marker='s')
        ps.append(p)

    plt.ylim(bottom=0, top=16)
    plt.grid(True)
    plt.legend(fontsize=7, loc='upper right', ncol=1)
    
    plt.subplot(gs[1, 0])

    # -- Expected error scaling of the two methods
    
    alpha_bse = -1
    alpha_dbse = -3

    alpha = -1 # scaling used in plot
    n_bse_1 = 3 # number of points to use for extrapolation fit
    n_bse_2 = 10 # number of points to use for extrapolation fit
    n_dbse = 3 # number of points to use for extrapolation fit

    # -- Fit DBSE
    
    nw = np.array([ float(p.nwf) for p in psd ])
    cG = np.array([ p.SzSz.chi_X.real for p in psd ])

    x = nw
    y = cG

    sidx = np.argsort(x)
    x, y = x[sidx], y[sidx]
        
    poly = np.polyfit(x[-n_dbse:]**alpha_dbse, y[-n_dbse:], 1)
    z_f = np.linspace(0, x[3]**alpha_dbse, num=1000)
    y_f = np.polyval(poly, z_f)
    
    plt.plot(z_f**(alpha/alpha_dbse), y_f, '-', color='gray', alpha=0.75, lw=1.0)

    for idx, X in enumerate(x):
        if X in skip_list:
            plt.plot(X**alpha, y[idx], 'o', alpha=0.75, color='gray', markersize=4)
        else:
            plt.plot(X**alpha, y[idx], 'o', alpha=0.75, color=colors[X], zorder=100, markersize=4)
    
    plt.plot([], [], 'o-', color='gray', label='DMFT+DBSE', alpha=0.75, lw=1.0, markersize=4)

    # -- Fit BSE

    nw = np.array([ float(p.nwf) for p in ps ])
    cG = np.array([ p.SzSz.chi_X.real for p in ps ])

    x = nw
    y = cG

    sidx = np.argsort(x)
    x, y = x[sidx], y[sidx]

    poly = np.polyfit(x[-n_bse_1:]**alpha_bse, y[-n_bse_1:], 1)
    z_f = np.linspace(0, x[3]**alpha_bse, num=1000)
    y_f = np.polyval(poly, z_f)
    plt.plot(z_f**(alpha/alpha_bse), y_f, '--', color='gray', alpha=0.75, lw=1.0)

    poly = np.polyfit(x[-n_bse_2:]**alpha_bse, y[-n_bse_2:], 2)
    z_f = np.linspace(0, x[3]**alpha_bse, num=1000)
    y_f = np.polyval(poly, z_f)
    plt.plot(z_f**(alpha/alpha_bse), y_f, ':', color='gray', alpha=0.75, lw=1.0)
    
    for idx, X in enumerate(x):
        if X in skip_list:
            plt.plot(X**alpha, y[idx], 's', alpha=0.75, color='gray', markersize=4)
        else:
            plt.plot(X**alpha, y[idx], 's', alpha=0.75, color=colors[X], zorder=100, markersize=4)

    
    plt.plot([], [], 's--', color='gray', label='DMFT+BSE', alpha=0.75, lw=1.0, markersize=4)

    plt.xlabel(r'$N_\nu^{' + f'{alpha}' + r'}$')
    plt.xlim(left=-0.001)
    plt.grid(True, axis='x')

    plt.xlim(left=-0.005)
    plt.ylim(bottom=0.0, top=4.5)
    plt.ylabel(r'$\chi_{S_z S_z}(\mathbf{q}_X)$')
    
    plt.savefig('figure_sro_chi_bandpath.svg')
    plt.show()
