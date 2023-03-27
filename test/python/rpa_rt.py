
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

from triqs_tprf.rpa_rt import chi0_tr_from_g_tr_les_gtr
from triqs_tprf.rpa_rt import chi0_tr_les_gtr_from_g_tr_les_gtr

from triqs_tprf.lattice import dynamical_screened_interaction_W

# Gregory coefficients of order 32

greg = np.array([
-7.8024383965700760956707426885259337723255157470703125E-01,
+1.7203004564829256217706188181182369589805603027343750E+00,
-1.1293914328292293802746826258953660726547241210937500E+01,
+6.9040528844601595892527257092297077178955078125000000E+01,
-3.5287838015253959156325436197221279144287109375000000E+02,
+1.5001985068704459536093054339289665222167968750000000E+03,
-5.3601248928468357917154207825660705566406250000000000E+03,
+1.6285104073565818907809443771839141845703125000000000E+04,
-4.2508001142523455200716853141784667968750000000000000E+04,
+9.6140272797904486651532351970672607421875000000000000E+04,
-1.8969528355236799689009785652160644531250000000000000E+05,
+3.2828785580717795528471469879150390625000000000000000E+05,
-5.0036565636323939543217420578002929687500000000000000E+05,
+6.7370192565594415646046400070190429687500000000000000E+05,
-8.0294330893410148564726114273071289062500000000000000E+05,
+8.4807946549424529075622558593750000000000000000000000E+05,
-7.9404676912516623269766569137573242187500000000000000E+05,
+6.5866466407562396489083766937255859375000000000000000E+05,
-4.8334433653750858502462506294250488281250000000000000E+05,
+3.1302010741073323879390954971313476562500000000000000E+05,
-1.7827630392221597139723598957061767578125000000000000E+05,
+8.8870604585172768565826117992401123046875000000000000E+04,
-3.8533790420769262709654867649078369140625000000000000E+04,
+1.4414059713544615078717470169067382812500000000000000E+04,
-4.6019186226658493978902697563171386718750000000000000E+03,
+1.2363601752227505130576901137828826904296875000000000E+03,
-2.7421022400211342073816922493278980255126953125000000E+02,
+4.8877487889135068144241813570261001586914062500000000E+01,
-6.7303812118926353846859456098172813653945922851562500E+00,
+6.7198303653292179582479093369329348206520080566406250E-01,
-4.3290258602310616176112034736433997750282287597656250E-02,
+1.3509659123127626401128509314730763435363769531250000E-03,
    ])

def fourier_trasform_of_spectral_function(t, f, g_f):

    Nt = len(t)
    Nf = len(f)
    zero_padding = (2*Nf) // Nt - 1
    df = f[1] - f[0]

    print(f'Nt = {len(t)}, Nf = {len(f)}, zero_padding = {zero_padding}')
    
    A_f = - g_f.imag / np.pi

    norm_trapz = np.trapz(A_f, x=f)
    print(f'norm_trapz = {norm_trapz}')

    norm_riemann = np.sum(A_f) * df
    print(f'norm_riemann = {norm_riemann}')

    #N = 2 * (2 - zero_padding) * Nf
    N = 4 * Nf
    print(f'N = {N}')
    window = np.blackman(Nf-1)
    window[:] = 1.
    #window = 1.

    norm = np.trapz(A_f[1:] * window, x=f[1:])
    print(f'norm = {norm}')

    #A_f /= norm

    
    #A_t = -1.j * np.fft.ifft(np.fft.ifftshift(A_f[:-1] * window, axes=0), n=N, axis=0)
    #A_t = -1.j * np.fft.fft(A_f[1:] * window, n=N, axis=0)
    A_t = -1.j * np.fft.fft(A_f[1:] * window, n=N, axis=0)

    N2 = N // (2 * (1 + zero_padding))
    A_t = A_t[:N2]

    #A_t /= (2. * np.pi)**2
    A_t *= df

    # -- Boundary corrections using Gregory weights

    #f = df * np.arange(
    
    k = len(greg)
    mat_tf = np.exp(1j*f[None, 1:k+1]*t[:, None]) * df * greg[None, :]
    #A_t += np.tensordot(mat_tf, A_f[1:k+1], axes=(1, 0))

    mat_tf = np.exp(1j*f[None, -k-1:-1]*t[:, None]) * df * greg[None, :]
    #A_t += np.tensordot(mat_tf, A_f[-k-1:-1], axes=(1, 0))
    
    t_shift = np.exp(1.j * t * f[-1])
    A_t *= t_shift
    
    if False:
        import matplotlib.pyplot as plt
        subp = [2, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        #plt.plot(t, t_shift.imag)
        plt.plot(t, A_t.real)
        plt.plot(t, A_t.imag)
        #plt.xlim([0, 10])
        plt.subplot(*subp); subp[-1] += 1
        #plt.plot(f, A_f * window)
        plt.plot(A_f)
        plt.plot(np.fft.fftshift(A_f[:-1]))
        plt.show(); exit()
    
    return A_t
    
def fourier_from_fX_to_tX(G_fX):
    pass

def test_fft(verbose=True):

    #for zero_padding in [0, 1, 2, 3]:
    for zero_padding in [0]:

        print(f'zero_padding = {zero_padding}')

        #dt = 0.0002
        #Nt = 2**16

        dt = 0.04 / 2**10
        Nt = 2**10 * 2**10

        t = dt * np.arange(Nt)
        tmesh = MeshReTime(t[0], t[-1], len(t))
        fmesh = fmesh_from_tmesh(tmesh, zero_padding)

        t = np.array(list(tmesh.values()))
        f = np.array(list(fmesh.values()))

        e = 0 - 0.5j
        #e = 1. - 0.5j
        g_t = -1j * np.exp(-1j * e * t)

        g_f = fourier_transformation_with_heaviside_matrix(f, t, g_t)
        g_f_ref = 1/(f - e)

        diff_f = np.max(np.abs(g_f - g_f_ref))
        print(f'Nt = {len(t)}, Nf = {len(f)}, zero_padding = {zero_padding}')
        print(f'diff_f = {diff_f:2.2E}')
        #np.testing.assert_array_almost_equal(g_f, g_f_ref, decimal=7)

        G_t = Gf(mesh=tmesh, target_shape=[1])
        G_t.data[:, 0] = -1j * np.exp(-1j * e * t)
        G_f = fourier_from_tX_to_fX(
            G_t, zero_padding=zero_padding, windowing=False)

        assert( fmesh == G_f.mesh )

        diff = np.max(np.abs(np.squeeze(G_f.data) - g_f_ref))
        print(f'diff = {diff:2.2E}')
        #np.testing.assert_array_almost_equal(
        #    np.squeeze(G_f.data), g_f_ref, decimal=7)

        # -- Frequency to time
        g_t_ref = fourier_trasform_of_spectral_function(t, f, g_f)
        #G_t_ref = fourier_from_fX_to_tX(G_f)

        diff_t = np.max(np.abs(g_t - g_t_ref))
        print(f'diff_t = {diff_t:2.2E}')
        
        

    if verbose:
        
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 9))

        subp = [5, 2, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, g_t.real)
        plt.plot(t, g_t_ref.real, '--')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, g_t.imag)
        plt.plot(t, g_t_ref.imag, '--')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, g_t.real - g_t_ref.real)

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(t, g_t.imag - g_t_ref.imag)
        
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


def test_chi0(verbose=True):

    beta = 100.0
    mu = 0.0
    nk = 8

    norb = 1
    a0 = 1.0
    kappa = 20.0
    t = 1.0

    Nt = 256 * 2
    dt = 0.1 / 2
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

    #del g_tk_les
    #del g_tk_gtr
    
    tmr = time.time()
    print('--> chi0_tr from g')
    chi0_tr = chi0_Tr_from_g_Tr_PH(g_tr_les, g_tr_gtr)
    print(f'done {time.time() - tmr} s')

    chi0_tr_ref = chi0_tr_from_g_tr_les_gtr(g_tr_les, g_tr_gtr)
    np.testing.assert_array_almost_equal(chi0_tr.data, chi0_tr_ref.data)

    chi0_tr_les, chi0_tr_gtr = chi0_tr_les_gtr_from_g_tr_les_gtr(g_tr_les, g_tr_gtr)
    chi0_tr_ref2 = chi0_tr_les - chi0_tr_gtr
    np.testing.assert_array_almost_equal(chi0_tr.data, chi0_tr_ref2.data)
    
    #del g_tr_les
    #del g_tr_gtr
    
    tmr = time.time()
    print('--> fourier_Tr_to_Tk (chi0)')
    chi0_tk = fourier_Tr_to_Tk(chi0_tr)
    print(f'done {time.time() - tmr} s')

    chi0_tk_les = fourier_Tr_to_Tk(chi0_tr_les)
    chi0_tk_gtr = fourier_Tr_to_Tk(chi0_tr_gtr)
    
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

    pi_tk_les = chi0_tk_les.copy()
    pi_tk_les.data[:] = -1.0 * pi_tk_les.data[:]

    pi_tk_gtr = chi0_tk_gtr.copy()
    pi_tk_gtr.data[:] = -1.0 * pi_tk_gtr.data[:]
    
    pi_fk = chi0_fk.copy()
    pi_fk.data[:] = -1.0 * pi_fk.data[:]

    print(f'--> W_fk'); tmr = time.time()
    W_fk_ref = dynamical_screened_interaction_W(pi_fk, V_k)
    print(f'done {time.time() - tmr} s')    

    tmr = time.time()
    print('--> RPA W_tk')
    Wc_tk = W_tk_RPA(pi_tk, V_k)
    print(f'done {time.time() - tmr} s')

    Wc_tk_les = W_tk_RPA(pi_tk_les, V_k)
    Wc_tk_gtr = W_tk_RPA(pi_tk_gtr, V_k)

    Wc_tk_ref = Wc_tk_gtr - Wc_tk_les
    np.testing.assert_array_almost_equal(Wc_tk.data, Wc_tk_ref.data)
    exit()
    
    tmr = time.time()
    print('--> fourier_from_tk_to_tr Wc')
    Wc_tr = fourier_Tk_to_Tr(Wc_tk)
    print(f'done {time.time() - tmr} s')

    Wc_tr_les = fourier_Tk_to_Tr(Wc_tk_les)
    Wc_tr_gtr = fourier_Tk_to_Tr(Wc_tk_gtr)
    
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

    from triqs_tprf.lattice import g0w_sigma
    print(f'--> sigma_fk'); tmr = time.time()
    sigma_fk_ref = g0w_sigma(mu, beta, e_k, W_fk, V_k, delta)
    print(f'done {time.time() - tmr} s')

    sigma_f_ref = interpolate_g_fk(sigma_fk_ref, k0)

    # -- Real time Sigma GW
    sigma_tr_les = g_tr_les.copy()
    sigma_tr_les.data[:] *= 1j * Wc_tr_les.data[..., 0, 0]

    sigma_tr_gtr = g_tr_gtr.copy()
    sigma_tr_gtr.data[:] *= 1j * Wc_tr_gtr.data[..., 0, 0]
    
    sigma_tr = sigma_tr_gtr - sigma_tr_les

    sigma_tk = fourier_Tr_to_Tk(sigma_tr)
    sigma_fk = fourier_from_tX_to_fX(
            sigma_tk, zero_padding=zero_padding, windowing=False)

    #sigma_fk.data[:] += V_k.data[None, :, :, :, 0, 0] * g_tk_les.data[0, ...][None, ...]
    #sigma_fk.data[:] *= 0.
    
    sigma_f = interpolate_g_fk(sigma_fk, k0)
    
    if verbose:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 8))
        
        subp = [4, 1, 1]

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

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(f, sigma_f.real, label='Re')
        plt.plot(f, sigma_f.imag, label='Im')
        plt.plot(f, sigma_f_ref.real, '--', label='Re', alpha=0.75)
        plt.plot(f, sigma_f_ref.imag, '--', label='Im', alpha=0.75)
        plt.ylabel(r'$\Sigma(k_{small}, \omega)$')
        plt.legend(loc='best')

        plt.tight_layout()
        plt.show()

    
if __name__ == "__main__":
        
    #test_fft()
    test_chi0()
