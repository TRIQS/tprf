"""

Author: Hugo U. R. Strand (2022)
"""


import time
import numpy as np

from triqs.gf import Gf, MeshProduct, MeshReFreq, MeshReTime, Idx

from triqs_tprf.lattice import fourier_fr_to_fk
from triqs_tprf.lattice import fourier_fk_to_fr
from triqs_tprf.lattice import fourier_Tr_to_Tk
from triqs_tprf.lattice import fourier_Tk_to_Tr
from triqs_tprf.lattice import g0_Tk_les_gtr_from_e_k
from triqs_tprf.lattice import chi0_Tr_from_g_Tr_PH


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


def fourier_from_tX_to_fX(g_tX, zero_padding=0, windowing=True):
    """ Fourier transform from real time to real frequency
    of Triqs Greens function with a real time mesh as first mesh.

    Accounts for Heavside function and gives both real and
    imaginary part of the real frequency Green's function."""

    in_mesh = g_tX.mesh
    is_product_mesh = type(in_mesh) == MeshProduct
    
    tmesh = in_mesh[0] if is_product_mesh else in_mesh
    assert( type(tmesh) == MeshReTime )
    
    fmesh = fmesh_from_tmesh(tmesh, zero_padding)
    
    if is_product_mesh:
        l = [fmesh] + [ in_mesh[i] for i in range(1, in_mesh.rank) ]
        out_mesh = MeshProduct(*l)
    else:
        out_mesh = fmesh
        
    g_fX = Gf(mesh=out_mesh, target_shape=g_tX.target_shape)
    
    t = np.array(list(tmesh.values()))
    f = np.array(list(fmesh.values()))

    g_fX.data[:] = fourier_transformation_with_heaviside_matrix(
        f, t, g_tX.data, windowing=windowing)

    return g_fX


def fmesh_from_tmesh(tmesh, zero_padding=0):
    """ Construct a real frequency mesh that is compatible
    with the given real time mesh (and the Fourier transforms). """
    
    assert(zero_padding >= 0)
    
    t = np.array(list(tmesh.values()))
    dt = t[1] - t[0]
    Nt = len(t)

    Nf = Nt * (1 + zero_padding) // 2
    f = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(Nf, d=dt*4))

    fmax = -f[0]
    df = f[1] - f[0]
    fmesh = MeshReFreq(omega_min=-fmax, omega_max=fmax-df, n_max=Nf)

    f_ref = np.array(list(fmesh))
    np.testing.assert_array_almost_equal(f, f_ref)
    
    return fmesh


def fourier_transformation_with_heaviside_matrix(
        w, t, G_t, windowing=False):

    """ Fourier transform G(t) to frequency assuming that G(t) = 0 for t < 0 

    G(w) = \int_0^\infty e^{i \omega t} G(t)

    NB! Making assumptions on the frequency and time grids.
    Use only with specifically constructed meshes. """
    
    Nt = len(t)
    Nf = len(w)
    zero_padding = (2*Nf) // Nt - 1

    #print(f'Nt = {len(t)}, Nw = {len(w)}, zero_padding = {zero_padding}')
    
    assert( Nt == G_t.shape[0] )
    
    if windowing:
        window = np.blackman(2*Nt)[Nt:].reshape([Nt] + [1]*len(G_t.shape[1:]))
    else:
        window = 1.
        
    # -- Fourier transform adding zeros for t < 0
    dt = t[1] - t[0]
    G_t = window * G_t
    N = 2*Nt*(1 + zero_padding)
    G_w = np.fft.ifft(G_t, n=N, axis=0) * N
    n = N // 8
    G_w = dt * np.fft.ifftshift(np.concatenate((G_w[:n], G_w[-n:])), axes=0)

    # -- Boundary corrections using Gregory weights
    k = len(greg)
    mat_wt = np.exp(1j*w[:, None]*t[None, :k]) * dt * greg[None, :]
    G_w += np.tensordot(mat_wt, G_t[:k], axes=(1, 0))
    
    return G_w


def chi0_tr_from_g_tr_les_gtr(g_tr_les, g_tr_gtr):

    rmesh = g_tr_les.mesh[-1]
    target_shape = list(g_tr_les.target_shape) * 2
    chi0_tr = Gf(mesh=g_tr_les.mesh, target_shape=target_shape)
    
    # Determine map of -r

    def minus_r_value(value, rmesh):
        value = -value
        for i, d in enumerate(rmesh.dims):
            value[i] = np.mod(value[i], d)
        return tuple(value)

    r_value_vec = [ tuple(r.value) for r in rmesh ]
    r_value_m_vec = [ minus_r_value(r.value, rmesh) for r in rmesh]
    r_m_idx = [ r_value_vec.index(r_m) for r_m in r_value_m_vec ]
    rmesh_list = list(rmesh)

    for r_p in rmesh:
        r_m = rmesh_list[r_m_idx[r_p.linear_index]] # look up mesh point corresponding to -r
        chi0_tr[:, r_p].data[:] = \
           1j * np.einsum('tda,tbc->tabcd', g_tr_les[:, r_p].data, np.conj(g_tr_gtr[:, r_m].data)) - \
           1j * np.einsum('tda,tbc->tabcd', g_tr_gtr[:, r_p].data, np.conj(g_tr_les[:, r_m].data))

    return chi0_tr


def chi0_tr_les_gtr_from_g_tr_les_gtr(g_tr_les, g_tr_gtr):

    mesh = g_tr_les.mesh
    rmesh = mesh[-1]
    target_shape = list(g_tr_les.target_shape) * 2
    chi0_tr_les = Gf(mesh=mesh, target_shape=target_shape)
    chi0_tr_gtr = Gf(mesh=mesh, target_shape=target_shape)
    
    # Determine map of -r

    def minus_r_value(value, rmesh):
        value = -value
        for i, d in enumerate(rmesh.dims):
            value[i] = np.mod(value[i], d)
        return tuple(value)

    r_value_vec = [ tuple(r.value) for r in rmesh ]
    r_value_m_vec = [ minus_r_value(r.value, rmesh) for r in rmesh]
    r_m_idx = [ r_value_vec.index(r_m) for r_m in r_value_m_vec ]
    rmesh_list = list(rmesh)

    for r_p in rmesh:

        # look up mesh point corresponding to -r
        r_m = rmesh_list[r_m_idx[r_p.linear_index]] 

        chi0_tr_les[:, r_p].data[:] = 1j * np.einsum(
            'tda,tbc->tabcd', g_tr_les[:, r_p].data, np.conj(g_tr_gtr[:, r_m].data))

        chi0_tr_gtr[:, r_p].data[:] = 1j * np.einsum(
            'tda,tbc->tabcd', g_tr_gtr[:, r_p].data, np.conj(g_tr_les[:, r_m].data))

    return chi0_tr_les, chi0_tr_gtr


def chi0_tk_from_ek(e_k, beta, tmesh):

    tmr = time.time()
    print('--> g_tk_les, g_tk_gtr')
    g_tk_les, g_tk_gtr = g0_Tk_les_gtr_from_e_k(e_k, tmesh, beta)
    print(f'done {time.time() - tmr} s')

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

    del chi0_tr
    
    return chi0_tk

    
def chi0_fk_from_ek(e_k, beta, tmesh, zero_padding=0, windowing=True):
    
    chi0_tk = chi0_tk_from_ek(e_k, beta, tmesh)

    tmr = time.time()
    print('--> fourier_from_tX_to_fX (chi0)')
    chi0_fk = fourier_from_tX_to_fX(
            chi0_tk, zero_padding=zero_padding, windowing=windowing)
    print(f'done {time.time() - tmr} s')
    
    return chi0_fk


def interpolate_g_k(g_k, k_vecs):
    assert( k_vecs.shape[1] == 3 )
    g_interp_k = np.zeros(
        [k_vecs.shape[0]] + list(g_k.target_shape), dtype=complex)

    for kidx, (kx, ky, kz) in enumerate(k_vecs):
        g_interp_k[kidx] = g_k(tuple([kx, ky, kz]))

    return np.squeeze(g_interp_k)


def interpolate_g_fk(g_fk, k_vecs):

    fmesh = g_fk.mesh[0]
    g_interp_fk = np.zeros((len(fmesh), k_vecs.shape[0]), dtype=complex)

    for fidx, f in enumerate(fmesh):
        g_k = g_fk[tuple([f, slice(None)])]
        g_interp_fk[fidx, :] = interpolate_g_k(g_k, k_vecs)

    return g_interp_fk


def W_tk_RPA_numpy(t, pi_tk, v_k):
    """ TODO: Implement PH rank 4 tensor products """
    
    W_tk = np.zeros_like(pi_tk)

    Nt = pi_tk.shape[0]
    dt = t[1] - t[0]

    Nk = pi_tk.shape[1]
    v_K = v_k.reshape([1] + list(v_k.shape))

    t0 = 0
    W_tk[t0, :] = v_K * pi_tk[t0, :] * v_K

    denom = 1. / (1. + dt*0.5*pi_tk[t0, :] * v_K)
    
    t1  = 1
    v_pi_v = v_K * pi_tk[t1, :].data * v_K
    I_K = dt * 0.5 * pi_tk[t1,:] *v_K * W_tk[t0,:]
    W_tk[t1, :] = (v_pi_v + I_K) * denom

    import tqdm
    for idx in tqdm.tqdm(range(2, Nt)):
        
        # -- convolution
        W_TK = W_tk[:idx, :]
        pi_TK = pi_tk[1:idx+1, :]
        
        I_TK = pi_TK[::-1] * v_K * W_TK
        I_K = np.sum(I_TK, axis=0) * dt
        
        I_K_boundary_corr = - 0.5*dt * I_TK[0]
        I_K += I_K_boundary_corr
            
        # -- RPA equation solution
        v_pi_v = v_K * pi_tk[idx, :] * v_K
        W_tk[idx, :] = (v_pi_v + I_K) * denom

    return W_tk


def W_tk_RPA(pi_tk, v_k):

    # -- Only supporting scalar interaction currently
    assert( v_k.target_shape == (1, 1, 1, 1) )
    
    tmesh, kmesh = pi_tk.mesh[0], pi_tk.mesh[1]
    t = np.array(list(tmesh.values()))

    V_k = v_k.data.reshape((len(kmesh), 1, 1, 1, 1))
    
    W_tk = pi_tk.copy()    
    W_tk.data[:] = W_tk_RPA_numpy(t, pi_tk.data, V_k)
    return W_tk
