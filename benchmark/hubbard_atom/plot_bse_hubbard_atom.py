# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import *

# ----------------------------------------------------------------------

from triqs_tprf.hubbard_atom import chi_ph_magnetic
from triqs_tprf.hubbard_atom import gamma_ph_magnetic
from triqs_tprf.hubbard_atom import single_particle_greens_function

# ----------------------------------------------------------------------

from triqs_tprf.linalg import inverse_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.chi_from_gg2 import chi_from_gg2_PH
from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued

# ----------------------------------------------------------------------
class Dummy():
     def __init__(self):
        pass

# ----------------------------------------------------------------------
def analytic_solution(beta, U, nw, nwf):

    data = Dummy()

    g_iw = single_particle_greens_function(beta=beta, U=U, nw=np.max([nw, nwf])*4)
    data.G_iw = g_iw

    # make block gf of the single gf
    G_iw_block = BlockGf(name_list=['up', 'dn'], block_list=[g_iw, g_iw])
    g_mat = block_iw_AB_to_matrix_valued(G_iw_block)
    
    # CHANGE NAME of the function to g2_ph_m!
    data.g2_ph_m = chi_ph_magnetic(beta=beta, U=U, nw=nw, nwf=nwf)

    data.chi_m = chi_from_gg2_PH(g_mat, data.g2_ph_m)
    data.chi0_m = chi0_from_gg2_PH(g_mat, data.g2_ph_m)    

    data.gamma_m = gamma_ph_magnetic(beta=beta, U=U, nw=nw, nwf=nwf)

    # -- Numeric BSE
    data.inv_chi0 = inverse_PH(data.chi0_m)
    #data.inv_chi = inverse_PH(data.chi_m)
    data.inv_chi = inverse_PH(data.g2_ph_m) # DEBUG TEST

    data.gamma_m_num = data.inv_chi0 - data.inv_chi
    
    data.label = r'Analytic'
    return data

# ----------------------------------------------------------------------
def fixed_fermionic_window(g2, nwf):

    nw = (g2.data.shape[0] + 1) / 2
    beta = g2.mesh.components[0].beta

    mesh_iw = MeshImFreq(beta=beta, S='Boson', n_max=nw)
    mesh_inu = MeshImFreq(beta=beta, S='Fermion', n_max=nwf)
    mesh_prod = MeshProduct(mesh_iw, mesh_inu, mesh_inu)

    g2_out = Gf(mesh=mesh_prod, target_shape=[1,1,1,1])

    n = g2.data.shape[1]
    s = n/2 - nwf
    e = n/2 + nwf
    
    g2_out.data[:] = g2.data[:, s:e, s:e]
    return g2_out

# ----------------------------------------------------------------------
if __name__ == '__main__':

    ana = analytic_solution(beta=2.0, U=5.0, nw=10, nwf=150)

    for key in ['gamma_m_num', 'gamma_m']:

        val = getattr(ana, key)
        val_fix = fixed_fermionic_window(val, nwf=10)
        setattr(ana, key, val_fix)

    print ana.gamma_m.data.shape
        
    for data in [ana.gamma_m.data, ana.gamma_m_num.data]:
        print np.max(np.abs(data.imag))
        np.testing.assert_array_almost_equal(data.imag, np.zeros_like(data))
    
    idx_list = [
        (0, 0, 0),
        (7, 10, 10),
        (9, 9, 9),
    #    (5, 10, 10),
    #    (5, 30, 30),
        ]
    for idx in idx_list:
        i1, i2, i3 = idx
        val = ana.gamma_m_num.data[i1, i2, i3]
        val_ref = ana.gamma_m.data[i1, i2, i3]
        print idx, val, val_ref, val - val_ref

    abs_diff = np.abs( ana.gamma_m_num.data - ana.gamma_m.data )
    diff = np.max(abs_diff)
    print 'max(abs(diff)) =', diff

    idx = np.argwhere(np.abs(abs_diff - diff) < 1e-9)
    print idx
    #exit()
    
    diff = np.max(np.abs( ana.gamma_m_num.data.imag - ana.gamma_m.data.imag ))
    print 'max(abs(Im[diff])) =', diff
        
    #exit()
    
    from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt
    plt.figure(figsize=(12, 8))

    subp = [3, 6, 1]

    #idx = 20
    idx = 10

    #opt = dict(vmin=-10., vmax=10., cmap='PuOr')
    opt = dict(cmap='PuOr')
    
    if True:
        g2 = ana.gamma_m
        s = np.squeeze(g2.data[idx, :, :])
        label = 'gamma ana'
        
        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [i,:,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, idx, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,i,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [:,i,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, :, idx])

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,:,i]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [:,:,i]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
    if True:
        g2 = ana.gamma_m_num
        s = np.squeeze(g2.data[idx, :, :])
        label = 'gamma num'
        
        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [i,:,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, idx, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,i,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [:,i,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, :, idx])

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,:,i]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [:,:,i]')
        plt.pcolormesh(s.imag)
        plt.colorbar()

    #opt = dict(vmin=-1., vmax=1., cmap='PuOr')
    opt = dict(cmap='PuOr')
    if True:
        g1 = ana.gamma_m
        g2 = ana.gamma_m_num
        s1 = np.squeeze(g1.data[idx, :, :])
        s2 = np.squeeze(g2.data[idx, :, :])
        s = s1 - s2
        label = 'gamma diff'
        
        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [i,:,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s1 = np.squeeze(g1.data[:, idx, :])
        s2 = np.squeeze(g2.data[:, idx, :])
        s = s1 - s2

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,i,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [:,i,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s1 = np.squeeze(g1.data[:, :, idx])
        s2 = np.squeeze(g2.data[:, :, idx])
        s = s1 - s2

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,:,i]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + label + ' [:,:,i]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
    plt.tight_layout()
    plt.show()
