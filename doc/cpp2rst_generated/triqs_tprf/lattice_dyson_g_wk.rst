..
   Generated automatically by cpp2rst

.. highlight:: c


.. _lattice_dyson_g_wk:

lattice_dyson_g_wk
==================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_wk_t lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k,
   triqs_tprf::g_w_cvt sigma_w)       (1)

    triqs_tprf::g_wk_t lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k,
   triqs_tprf::g_wk_cvt sigma_wk)     (2)

(1)Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`


Parameters
----------

 * **mu**: chemical potential :math:`\mu`

 * **e_k**: discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

 * **sigma_w**: imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`



Return value
------------

Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$

Documentation
-------------


 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, and a momentum independent Matsubara frequency
 self energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`.



(2)Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`


Parameters
----------

 * **mu**: chemical potential :math:`\mu`

 * **e_k**: discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

 * **sigma_wk**: imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})`



Return value
------------

Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$

Documentation
-------------


 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, and a momentum independent Matsubara frequency
 self energy :math:`\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})`.