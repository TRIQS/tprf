
import time
import numpy as np

from triqs.gf import Gf, MeshProduct, MeshReFreq, MeshReTime, Idx
from triqs.lattice.tight_binding import TBLattice


from triqs_tprf.lattice import lattice_dyson_g0_fk
from triqs_tprf.lattice import lindhard_chi00


from triqs_tprf.rpa_rt import fmesh_from_tmesh
from triqs_tprf.rpa_rt import fourier_from_tX_to_fX
from triqs_tprf.rpa_rt import fourier_transformation_with_heaviside_matrix

#from triqs_tprf.rpa_rt import chi0_fk_from_ek
from triqs_tprf.rpa_rt import interpolate_g_k, interpolate_g_fk
from triqs_tprf.rpa_rt import W_tk_RPA

from triqs_tprf.lattice import g0_Tk_les_gtr_from_e_k
from triqs_tprf.lattice import fourier_Tr_to_Tk
from triqs_tprf.lattice import fourier_Tk_to_Tr
from triqs_tprf.lattice import chi0_Tr_from_g_Tr_PH

from triqs_tprf.lattice import dynamical_screened_interaction_W


def test_fft(verbose=False):

    for zero_padding in [0, 1, 2, 3]:

        print(f'zero_padding = {zero_padding}')

        dt = 0.04
        Nt = 2**10

        t = dt * np.arange(Nt)
        tmesh = MeshReTime(t[0], t[-1], len(t))
        fmesh = fmesh_from_tmesh(tmesh, zero_padding)

        t = np.array(list(tmesh.values()))
        f = np.array(list(fmesh.values()))

        e = 1 - 0.5j
        g_t = -1j * np.exp(-1j * e * t)

        g_f = fourier_transformation_with_heaviside_matrix(f, t, g_t)
        g_f_ref = 1/(f - e)

        diff = np.max(np.abs(g_f - g_f_ref))
        print(f'Nt = {len(t)}, Nf = {len(f)}, zero_padding = {zero_padding}')
        print(f'diff = {diff:2.2E}')
        np.testing.assert_array_almost_equal(g_f, g_f_ref, decimal=7)

        G_t = Gf(mesh=tmesh, target_shape=[1])
        G_t.data[:, 0] = -1j * np.exp(-1j * e * t)
        G_f = fourier_from_tX_to_fX(
            G_t, zero_padding=zero_padding, windowing=False)

        assert( fmesh == G_f.mesh )

        diff = np.max(np.abs(np.squeeze(G_f.data) - g_f_ref))
        print(f'diff = {diff:2.2E}')
        np.testing.assert_array_almost_equal(
            np.squeeze(G_f.data), g_f_ref, decimal=7)

    if verbose:
        
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 10))

        subp = [4, 2, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, g_t.real)

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, g_t.imag)

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, np.abs(g_t.real))
        plt.semilogy([], [])

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, np.abs(g_t.imag))
        plt.semilogy([], [])

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, g_f.real)
        plt.plot(f, g_f_ref.real, '--')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, g_f.imag)
        plt.plot(f, g_f_ref.imag, '--')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, np.abs(g_f.real - g_f_ref.real))
        plt.semilogy([], [])

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, np.abs(g_f.imag - g_f_ref.imag))
        plt.semilogy([], [])

        plt.tight_layout()
        plt.show()


def test_chi0(verbose=False):

    beta = 100.0
    mu = 0.0
    nk = 8

    norb = 1
    a0 = 1.0
    kappa = 20.0
    t = 1.0

    Nt = 256
    dt = 0.1
    zero_padding = 0
    delta = 0.4

    units = [(a0, 0, 0),
             (0, a0, 0)]

    hop= {(+1,0) : t * np.eye(norb),       
          (-1,0) : t * np.eye(norb),     
          (0,+1) : t * np.eye(norb),
          (0,-1) : t * np.eye(norb)}

    orb_pos = [(0,0,0)] * norb

    H = TBLattice(units, hop, orb_pos)

    kmesh = H.get_kmesh(n_k=(nk, nk, 1))
    e_k = H.fourier(kmesh)
    e_k.data[:] -= -2.0

    t = dt * np.arange(Nt)
    tmesh = MeshReTime(t[0], t[-1], len(t))
        
    tmr = time.time()
    print('--> g_tk_les, g_tk_gtr')
    g_tk_les, g_tk_gtr = g0_Tk_les_gtr_from_e_k(e_k, tmesh, beta)
    print(f'done {time.time() - tmr} s')

    g_tk_les.data[:] *= np.exp(- 0.5*delta * t)[:, None, None, None]
    g_tk_gtr.data[:] *= np.exp(- 0.5*delta * t)[:, None, None, None]

    tmr = time.time()
    print('--> fourier_from_Tk_to_Tr (g)')
    g_tr_les = fourier_Tk_to_Tr(g_tk_les)
    g_tr_gtr = fourier_Tk_to_Tr(g_tk_gtr)
    print(f'done {time.time() - tmr} s')

    del g_tk_les
    del g_tk_gtr
    
    tmr = time.time()
    print('--> chi0_tr from g')
    chi0_tr = chi0_Tr_from_g_Tr_PH(g_tr_les, g_tr_gtr)
    print(f'done {time.time() - tmr} s')

    del g_tr_les
    del g_tr_gtr
    
    tmr = time.time()
    print('--> fourier_Tr_to_Tk (chi0)')
    chi0_tk = fourier_Tr_to_Tk(chi0_tr)
    print(f'done {time.time() - tmr} s')

    tmr = time.time()
    print('--> fourier_from_tX_to_fX (chi0)')
    chi0_fk = fourier_from_tX_to_fX(
            chi0_tk, zero_padding=zero_padding, windowing=False)
    print(f'done {time.time() - tmr} s')
    
    fmesh = chi0_fk.mesh[0]
    
    print(f'--> lattice_dyson_g0_fk'); tmr = time.time()
    g0_fk = lattice_dyson_g0_fk(mu, e_k, fmesh, delta)
    print(f'done {time.time() - tmr} s')

    print(f'--> lindhard_chi00'); tmr = time.time()
    chi0_fk_ref = lindhard_chi00(e_k, fmesh, beta, mu, delta)
    print(f'done {time.time() - tmr} s')
    
    kx, ky, kz = [0.1]*3
    k0 = np.array([[kx, ky, kz]])

    chi0_f = interpolate_g_fk(chi0_fk, k0)
    chi0_f_ref = interpolate_g_fk(chi0_fk_ref, k0)
    f = np.array(list(fmesh.values()))

    chi0_diff = np.max(np.abs(chi0_f - chi0_f_ref))
    print(f'chi0_diff = {chi0_diff:2.2E}')
    np.testing.assert_array_almost_equal(chi0_f, chi0_f_ref, decimal=5)

    # -- Add W ?

    # -- Coulomb interaction

    def backfold(kvec, kmesh):
        # get the q vector in internal units
        kvecInt = kvec @ np.linalg.inv(kmesh.domain.units)

        # back-fold if nececcary
        for ii in range(len(kvecInt)):
            if(kvecInt[ii] >= 0.5):
                kvecInt[ii] -= 1.0 
            if(kvecInt[ii] < -0.5):
                kvecInt[ii] += 1.0 

        # get back-folded vector in 1/ang
        kvecFolded = kvecInt @ kmesh.domain.units
        return kvecFolded

    V_k = Gf(mesh=kmesh, target_shape=[1]*4)
    for k in kmesh:
        kind = k.linear_index
        kfolded = backfold(k.value, kmesh)
        knorm = np.linalg.norm(kfolded)

        if(np.isclose(knorm, 0.0)): 
            V_k.data[kind] = 0.0
            continue

        e2 = 14.399
        V_k.data[kind] = 2.0 * np.pi * e2 / (knorm * kappa * a0*a0)

    
    pi_tk = chi0_tk.copy()
    pi_tk.data[:] = -1.0 * pi_tk.data[:]

    pi_fk = chi0_fk.copy()
    pi_fk.data[:] = -1.0 * pi_fk.data[:]

    print(f'--> W_fk'); tmr = time.time()
    W_fk_ref = dynamical_screened_interaction_W(pi_fk, V_k)
    print(f'done {time.time() - tmr} s')    

    tmr = time.time()
    print('--> RPA W_tk')
    Wc_tk = W_tk_RPA(pi_tk, V_k)
    print(f'done {time.time() - tmr} s')

    tmr = time.time()
    print('--> fourier_from_tk_to_fk chi0')
    Wc_fk = fourier_from_tX_to_fX(Wc_tk, zero_padding=zero_padding, windowing=False)
    W_fk = Wc_fk.copy()
    W_fk.data[:] += V_k.data
    print(f'done {time.time() - tmr} s')

    W_f = interpolate_g_fk(W_fk, k0)
    W_f_ref = interpolate_g_fk(W_fk_ref, k0)

    W_diff = np.max(np.abs(W_f - W_f_ref))
    print(f'W_diff = {W_diff:2.2E}')
    np.testing.assert_array_almost_equal(W_f, W_f_ref, decimal=2)
    
    if verbose:
        import matplotlib.pyplot as plt
        subp = [3, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, chi0_f.real, label='Re')
        plt.plot(f, chi0_f.imag, label='Im')
        plt.plot(f, chi0_f_ref.real, '--', label='Re', alpha=0.75)
        plt.plot(f, chi0_f_ref.imag, '--', label='Im', alpha=0.75)
        plt.ylabel(r'$chi_0(k_{small}, \omega)$')
        plt.legend(loc='best')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, W_f.real, label='Re')
        plt.plot(f, W_f.imag, label='Im')
        plt.plot(f, W_f_ref.real, '--', label='Re', alpha=0.75)
        plt.plot(f, W_f_ref.imag, '--', label='Im', alpha=0.75)
        plt.ylabel(r'$W(k_{small}, \omega)$')
        plt.legend(loc='best')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, np.abs(W_f - W_f_ref))
        plt.ylabel(r'$\Delta W(k_{small}, \omega)$')
        plt.legend(loc='best')
        
        plt.show()

    
if __name__ == "__main__":
        
    test_fft(verbose=True)
    test_chi0(verbose=True)
