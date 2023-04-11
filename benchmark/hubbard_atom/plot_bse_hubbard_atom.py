# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import *

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.analytic_hubbard_atom import analytic_hubbard_atom

# ----------------------------------------------------------------------
def fixed_fermionic_window(g2, nwf):

     from cpp2py.compiler import compile

     code = """

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/utility/timer.hpp>

using namespace triqs;   
using namespace triqs::gfs;
using namespace triqs::mesh;

typedef gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_iw_t;

g2_iw_t fixed_fermionic_window(g2_iw_t g2_in, int nwf) {

     auto mesh_b = std::get<0>(g2_in.mesh());

     mesh::imfreq mesh_f{mesh_b.beta, Fermion, nwf};

     g2_iw_t g2_out{{mesh_b, mesh_f, mesh_f}, {1, 1, 1, 1}};

     clef::placeholder<0> w;
     clef::placeholder<1> n1;
     clef::placeholder<2> n2;

     g2_out(w, n1, n2) << g2_in(w, n1, n2);

     return g2_out;
}
     """
     
     cxxflags = '-O3 -march=native -Ofast -mavx -mfma -mavx2 -g -ggdb -Wno-invalid-partial-specialization '
     includes = ' -I /opt/local/include/ -I /opt/local/include/openmpi-clang50'
     M = compile(code, modules = "triqs", cxxflags=cxxflags + includes)

     g2_out = M.fixed_fermionic_window(g2, nwf)

     return g2_out

# ----------------------------------------------------------------------
def fixed_fermionic_window_python(g2, nwf):

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
def max_abs_diff_in_window(g2_a, g2_b, nwf):

    g2_a = fixed_fermionic_window(g2_a, nwf)
    g2_b = fixed_fermionic_window(g2_b, nwf)

    print(g2_a.data.shape)
    print(g2_b.data.shape)
    
    return np.max(np.abs( g2_a.data - g2_b.data ))
     
# ----------------------------------------------------------------------
def window_conv_depr():

    nwf_vec = np.array([5, 10, 20, 40, 80, 160, 320])
    diff_vec = np.zeros_like(nwf_vec, dtype=np.float)

    for idx, nwf in enumerate(nwf_vec):
         d = analytic_solution(beta=2.0, U=5.0, nw=1, nwf=nwf)
         diff = np.max(np.abs( d.gamma_m.data - d.gamma_m_num.data ))
         diff_vec[idx] = diff
         print('nwf, diff =', idx, nwf, diff)

    print(diff_vec)
         
    from triqs.plot.mpl_interface import oplot, oplotr, oploti, plt

    x = 1./nwf_vec

    plt.figure(figsize=(3.25, 3))
    
    plt.plot(x, diff_vec, 'o-', alpha=0.75)
    plt.xlabel(r'$1/n_{w}$')
    plt.ylabel(r'$\max|\Gamma_{ana} - \Gamma_{num}|$')
    plt.ylim([0, diff_vec.max()])
    plt.xlim([0, x.max()])

    plt.tight_layout()
    plt.savefig('figure_bse_hubbard_atom_convergence.pdf')
    plt.show()


# ----------------------------------------------------------------------
if __name__ == '__main__':

    p = ParameterCollection(beta=20., U=5., nw=1, nwf=20, nwf_gf=80)
    #ana = analytic_solution(**p.dict())
    ana = analytic_hubbard_atom(**p.dict())

    ana.gamma_m.data[:] -= p.U
    ana.gamma_m_num.data[:] -= p.U
    
    #for data in [ana.gamma_m.data, ana.gamma_m_num.data]:
    #    #print np.max(np.abs(data.imag))
    #    np.testing.assert_array_almost_equal(data.imag, np.zeros_like(data))
    
    diff = np.max(np.abs( ana.gamma_m_num.data.imag - ana.gamma_m.data.imag ))
    print('max(abs(Im[diff])) =', diff)

    diff = np.max(np.abs( ana.gamma_m_num.data.real - ana.gamma_m.data.real ))
    print('max(abs(Re[diff])) =', diff)
        
    #exit()
    
    from triqs.plot.mpl_interface import oplot, oplotr, oploti, plt
    plt.figure(figsize=(6*2, 8))

    #subp = [3, 6, 1]
    subp = [4, 1, 1]

    #idx = 20
    idx = 0

    vmax = np.max(np.abs(ana.gamma_m.data.real))    
    opt = dict(vmin=-vmax, vmax=vmax, cmap='PuOr')
    #opt = dict(vmin=-10., vmax=10., cmap='PuOr')
    #opt = dict(cmap='PuOr')

    ax = plt.subplot(*subp); subp[-1] += 1
    oplot(ana.G_iw)

    
    if True:
        g2 = ana.gamma_m
        label = 'gamma ana'

        s = np.squeeze(g2.data[0, :, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()        

    if False:
        g2 = ana.gamma_m
        label = 'gamma ana'

        s = np.squeeze(g2.data[idx, :, :])

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
        label = 'gamma num'

        s = np.squeeze(g2.data[0, :, :])

        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real, **opt)
        plt.colorbar()

    if False:
        g2 = ana.gamma_m_num
        label = 'gamma num'

        s = np.squeeze(g2.data[idx, :, :])

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
        s1 = np.squeeze(g1.data[0, :, :])
        s2 = np.squeeze(g2.data[0, :, :])
        s = s1 - s2
        label = 'gamma diff'

        vmax = np.max(np.abs(s.real))    
        opt_diff = dict(vmin=-vmax, vmax=vmax, cmap='PuOr')
        
        ax = plt.subplot(*subp); subp[-1] += 1
        plt.title('Re ' + label + ' [i,:,:]')
        plt.pcolormesh(s.real, **opt_diff)
        plt.colorbar()

    if False:
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
