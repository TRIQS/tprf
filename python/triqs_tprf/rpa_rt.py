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


def fermi_function(E, beta):
    """
    Evaluate the Fermi distribution function at energy :math:`E` and inverse temperature :math:`\beta`.

    The Fermi distribution function :math:`f_\\beta(E)` has the form

    .. math:: f_\\beta(E) = \\frac{1}{1 + e^{\\beta E}}

    the evaluation is stabilized using separate formulas for :math:`E \lessgtr 0`.

    Parameters
    ----------

    E : ndarray
        Energies

    beta : float
        Inverse temperature

    Returns
    -------

    f : ndarray
        The Fermi distribution function evaluated at :math:`E`
    
    """
    
    f = np.zeros_like(E)
    p, m = np.argwhere(E > 0), np.argwhere(E <= 0)

    f[p] = np.exp(-beta*E[p]) / (1. + np.exp(-beta*E[p]))
    f[m] = 1. / (np.exp(beta*E[m]) + 1.)

    return f


def fourier_transformation_with_heaviside_matrix(
        w, t, G_t, greg=None, order=16, windowing=False, kaiser_beta=14.):

    """ Fourier transform G(t) to frequency assuming that G(t) = 0 for t < 0 

    G(w) = \int_0^\infty e^{i \omega t} G(t)

    NB! Making assumptions on the frequency and time grids.
    Use only with specifically constructed meshes. """

    if greg is None:
        from pyneco.gregory import gregory_coeff    
        greg = -1 + np.array(gregory_coeff(order), dtype=float)

    n = G_t.shape[0]
    if windowing:
        #window = np.blackman(2*n)[n:].reshape([n] + [1]*len(G_t.shape[1:]))
        window = np.kaiser(2*n, kaiser_beta)[n:].reshape([n] + [1]*len(G_t.shape[1:]))
    else:
        window = 1.
        
    # -- Fourier transform adding zeros for t < 0
    dt = t[1] - t[0]
    G_t = window * G_t
    G_w = np.fft.ifft(G_t, n=n*2, axis=0) * G_t.shape[0]*2
    n = len(w) // 2
    G_w = dt * np.fft.ifftshift(np.concatenate((G_w[:n], G_w[-n:])), axes=0)

    # -- Boundary corrections using Gregory weights
    k = len(greg)
    mat_wt = np.exp(1j*w[:, None]*t[None, :k]) * dt * greg[None, :]
    G_w += np.tensordot(mat_wt, G_t[:k], axes=(1, 0))
    
    return G_w


def fourier_from_tk_to_fk(g_tk, g_fk, windowing=True, kaiser_beta=14.):
    """ Inplace FFT """
    tmesh = g_tk.mesh[0]
    fmesh = g_fk.mesh[0]
    t = np.array(list(tmesh.values()))
    f = np.array(list(fmesh.values()))
    g_fk.data[:] = fourier_transformation_with_heaviside_matrix(
        f, t, g_tk.data, windowing=windowing, kaiser_beta=kaiser_beta)


def setup_time_and_frequency_meshes(fmax, df, zero_padding):

    n = int(np.round(2*fmax / df))
    f = np.linspace(-fmax, fmax, num=n, endpoint=False)

    df = f[1] - f[0]
    Nf = len(f) * (1 + 2*zero_padding)
    dt = 2*np.pi/df/Nf
    Nt = Nf // 2
    t = dt * np.arange(Nt)

    tmesh = MeshReTime(t[0], t[-1], len(t))
    fmesh = MeshReFreq(omega_min=-fmax, omega_max=fmax-df, n_max=len(f))

    t_triqs = np.array(list(tmesh.values()))
    f_triqs = np.array(list(fmesh.values()))

    np.testing.assert_array_almost_equal(t, t_triqs)
    np.testing.assert_array_almost_equal(f, f_triqs)

    return tmesh, fmesh


def g_tk_ret_les_gtr(e_k, tmesh, beta):

    kmesh = e_k.mesh
    
    g_tk_les = Gf(mesh=MeshProduct(tmesh, kmesh), target_shape=e_k.target_shape)
    g_tk_gtr = Gf(mesh=MeshProduct(tmesh, kmesh), target_shape=e_k.target_shape)
    g_tk_ret = Gf(mesh=MeshProduct(tmesh, kmesh), target_shape=e_k.target_shape)

    for k in kmesh:
        H = e_k[k]
        E, U = np.linalg.eigh(H)
        H_ref = U @ np.diag(E) @ np.conj(U).T
        np.testing.assert_array_almost_equal(H_ref, H)

        f_E = fermi_function(E, beta)

        g_E_les = +1.j * f_E[None, :] * np.exp(-1j * E[None, :] * t[:, None])
        g_E_gtr = -1.j * (1 - f_E[None, :]) * np.exp(-1j * E[None, :] * t[:, None])

        g_tk_les[:, k].data[:] = np.einsum('ab,tb,bc->tac', U, g_E_les, np.conj(U).T) 
        g_tk_gtr[:, k].data[:] = np.einsum('ab,tb,bc->tac', U, g_E_gtr, np.conj(U).T)
        g_tk_ret[:, k].data[:] = g_tk_gtr[:, k].data[:] - g_tk_les[:, k].data[:]

    return g_tk_ret, g_tk_les, g_tk_gtr
    

def chi0_tr_from_g_tr_les_gtr(g_tr_les, g_tr_gtr):

    rmesh = g_tr_les.mesh[-1]
    target_shape = list(g_tr_les.target_shape) * 2
    chi0_tr = Gf(mesh=MeshProduct(tmesh, rmesh), target_shape=target_shape)
    
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


def chi0_from_ek(e_k, beta, tmesh, fmesh, kaiser_beta=14.):
    
    tmr = time.time()
    print('--> g_tk')
    g_tk_les, g_tk_gtr = g0_Tk_les_gtr_from_e_k(e_k, tmesh, beta)
    print(f'done {time.time() - tmr} s')

    tmr = time.time()
    print('--> fourier_from_Tk_to_Tr (g)')
    g_tr_les = fourier_Tk_to_Tr(g_tk_les)
    g_tr_gtr = fourier_Tk_to_Tr(g_tk_gtr)
    print(f'done {time.time() - tmr} s')
    
    tmr = time.time()
    print('--> chi0_tr from g')
    chi0_tr = chi0_Tr_from_g_Tr_PH(g_tr_les, g_tr_gtr)
    print(f'done {time.time() - tmr} s')
    
    tmr = time.time()
    print('--> fourier_Tr_to_Tk (chi0)')
    chi0_tk = fourier_Tr_to_Tk(chi0_tr)
    print(f'done {time.time() - tmr} s')
    
    # -- Transform to frequency

    tmr = time.time()
    print('--> fourier_from_tk_to_fk chi0')
    kmesh = chi0_tk.mesh[-1]
    chi0_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=chi0_tk.target_shape)
    fourier_from_tk_to_fk(chi0_tk, chi0_fk, kaiser_beta=kaiser_beta)
    print(f'done {time.time() - tmr} s')

    return chi0_fk
