# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import *
from pytriqs.operators import *

from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from triqs_tprf.hubbard_atom import chi_ph_magnetic
from triqs_tprf.hubbard_atom import gamma_ph_magnetic
from triqs_tprf.hubbard_atom import single_particle_greens_function

# ----------------------------------------------------------------------

from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued
from triqs_tprf.freq_conv import from_3nu_PH

from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH
from triqs_tprf.chi_from_gg2 import chi_from_gg2_PH

from triqs_tprf.linalg import inverse_PH
    
# ----------------------------------------------------------------------
class Dummy():
     def __init__(self):
        pass

# ----------------------------------------------------------------------
def read_pomerol_data(filename):

    data = Dummy()

    print '--> Loading:', filename
    with HDFArchive(filename, 'r') as h5:
        for key, value in h5.items():
            setattr(data, key, value)                

    a, b =  data.G2_iw_ph[('up', 'up')].data.shape[:2]
    nw, nwf = (a+1)/2, b/2
    data.nw, data.nwf = nw, nwf    
            
    data.beta = data.params['beta']
    data.U = data.params['U']

    # -- Single particle Green's function
            
    g_mat = block_iw_AB_to_matrix_valued(data.G_iw)
    g_mat.name = 'g_mat'

    data.G_iw = data.G_iw['up'] # For plotting just take one spin

    # -- Two particle Green's function

    # three ferimonic
    
    g2_uu = data.G2_iw_AABB[('up','up')]
    g2_ud = data.G2_iw_AABB[('up','dn')]

    g2_m = g2_uu - g2_ud
    
    # -- Two particle Green's function in PH channel

    beta = data.params['beta']
    mesh_iw = MeshImFreq(beta=beta, S='Boson', n_max=nw)
    mesh_inu = MeshImFreq(beta=beta, S='Fermion', n_max=nwf)
    mesh_prod = MeshProduct(mesh_iw, mesh_inu, mesh_inu)

    data.g2_ph_m = Gf(mesh=mesh_prod, target_shape=[1,1,1,1])

    from_3nu_PH(data.g2_ph_m, g2_m)

    # -- reference

    data.g2_ph_m_ref = data.G2_iw_ph[('up', 'up')] - data.G2_iw_ph[('up','dn')]
    
    # -- Generalized susceptibility in PH channel
    
    data.chi_m = chi_from_gg2_PH(g_mat, data.g2_ph_m)
    data.chi0_m = chi0_from_gg2_PH(g_mat, data.g2_ph_m)

    # -- Invert BSE

    data.gamma_m = inverse_PH(data.chi0_m) - inverse_PH(data.chi_m)

    data.label = r'Pomerol'

    return data
    
# ----------------------------------------------------------------------
def analytic_solution(beta, U, nw, nwf):

    data = Dummy()

    g_iw = single_particle_greens_function(beta=beta, U=U, nw=nw*2)
    data.G_iw = g_iw

    # make block gf of the single gf
    G_iw_block = BlockGf(name_list=['up', 'dn'], block_list=[g_iw, g_iw])
    g_mat = block_iw_AB_to_matrix_valued(G_iw_block)
    
    # CHANGE NAME of the function to g2_ph_m!
    data.g2_ph_m = chi_ph_magnetic(beta=beta, U=U, nw=nw, nwf=nwf)

    data.chi_m = chi_from_gg2_PH(g_mat, data.g2_ph_m)
    data.chi0_m = chi0_from_gg2_PH(g_mat, data.g2_ph_m)    

    data.gamma_m = gamma_ph_magnetic(beta=beta, U=U, nw=nw, nwf=nwf)
    data.gamma_m_num = inverse_PH(data.chi0_m) - inverse_PH(data.chi_m)
    
    data.label = r'Analytic'
    return data

# ----------------------------------------------------------------------
if __name__ == '__main__':

    filename = 'data_pomerol_hubbard_atom.h5'
    pom = read_pomerol_data(filename)
    np.testing.assert_array_almost_equal(pom.g2_ph_m.data, pom.g2_ph_m_ref.data)

    ana = analytic_solution(beta=pom.beta, U=pom.U, nw=pom.nw, nwf=pom.nwf)
    np.testing.assert_array_almost_equal(pom.G_iw.data, ana.G_iw.data)
    
    np.testing.assert_array_almost_equal(pom.g2_ph_m.data, ana.g2_ph_m.data)

    np.testing.assert_array_almost_equal(pom.chi_m.data, ana.chi_m.data)
    np.testing.assert_array_almost_equal(pom.chi0_m.data, ana.chi0_m.data)

    exit()
    
    for d in [pom, ana]:
        print d.label
        print d.G_iw.data.shape
        print d.chi_m.data.shape

        data = d.chi_m.data.imag
        #np.testing.assert_array_almost_equal(data, np.zeros_like(data)) 
        
    from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt

    if False:
        plt.figure()
        for d in [pom, ana]:
            oplotr(d.G_iw, label=d.label)
            oploti(d.G_iw, label=d.label)

    plt.figure(figsize=(12, 8))

    subp = [4, 6, 1]

    #idx = 20
    idx = 1

    print 'TEST ' * 10
    print 'replacing chi_m with g2_m'


    
    #pom.chi_m = pom.g2_ph_m_ref
    
    if True:
        g2 = pom.g2_ph_m_ref
        s = np.squeeze(g2.data[idx, :, :])
        label = 'g2_ph_m_ref'
        
        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, idx, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,i,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, :, idx])

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,:,i]')
        plt.pcolormesh(s.real)
        plt.colorbar()

    if True:
        g2 = pom.g2_ph_m
        s = np.squeeze(g2.data[idx, :, :])
        label = 'g2_ph_m'

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, idx, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,i,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()
        
        s = np.squeeze(g2.data[:, :, idx])

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [:,:,i]')
        plt.pcolormesh(s.real)
        plt.colorbar()
        
    for d in [pom, ana]:

        s = np.squeeze(d.chi_m.data[idx, :, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + d.label + ' [i,:,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('I, ' + d.label + ' [i,:,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s = np.squeeze(d.chi_m.data[:, idx, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + d.label + ' [:,i,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + d.label + ' [:,i,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s = np.squeeze(d.chi_m.data[:, :, idx])

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + d.label + ' [:,:,i]')
        plt.pcolormesh(s.real)
        plt.colorbar()

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + d.label + ' [:,:,i]')
        plt.pcolormesh(s.imag)
        plt.colorbar()


    if True:

        d1 = pom
        d2 = ana
        
        s1 = np.squeeze(d1.chi_m.data[idx, :, :])
        s2 = np.squeeze(d2.chi_m.data[idx, :, :])
        s = s1 - s2

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + d.label + ' [i,:,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('I, ' + d.label + ' [i,:,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s1 = np.squeeze(d1.chi_m.data[:, idx, :])
        s2 = np.squeeze(d2.chi_m.data[:, idx, :])
        s = s1 - s2

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + d.label + ' [:,i,:]')
        plt.pcolormesh(s.real)
        plt.colorbar()

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + d.label + ' [:,i,:]')
        plt.pcolormesh(s.imag)
        plt.colorbar()
        
        s1 = np.squeeze(d1.chi_m.data[:, :, idx])
        s2 = np.squeeze(d2.chi_m.data[:, :, idx])
        s = s1 - s2

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + d.label + ' [:,:,i]')
        plt.pcolormesh(s.real)
        plt.colorbar()

        a = plt.subplot(*subp); subp[-1] += 1
        plt.title('Im ' + d.label + ' [:,:,i]')
        plt.pcolormesh(s.imag)
        plt.colorbar()

        
    plt.tight_layout()
    plt.show()
