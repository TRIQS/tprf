..
   Generated automatically by cpp2rst

.. highlight:: c


.. _eliashberg_product_fft:

eliashberg_product_fft
======================

**Synopsis**:

.. code-block:: c

    triqs_tprf::gk_iw_t eliashberg_product_fft (triqs_tprf::chi_tr_vt Gamma_pp_dyn_tr,
   triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::gk_iw_vt g_wk,
   triqs_tprf::gk_iw_vt delta_wk)

Linearized Eliashberg product via FFT


Parameters
----------

 * **chi_rt**: dynamic part of the particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{r}, \tau)`

 * **chi_r**: constant part of the particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{r})`

 * **g_kw**: single particle Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\nu_n)`

 * **delta_kw**: pairing self-energy :math:`\Delta_{\bar{a}\bar{b}}(\mathbf{k}, i\nu_n)`



Return value
------------

Gives the result of the product :math:`\Delta^{(out)} \sim \Gamma^{(pp)}GG \Delta`

Documentation
-------------


     Computes the product

     .. math::
         \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{k}'} \sum_{i\nu'}
	 \Gamma_{A\bar{a}B\bar{b}}(\mathbf{k}-\mathbf{k}', i\nu - i\nu')
	 \\ \times
	 G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')\,,

     by taking advantage of the convolution theorem.

     We therefore first calculate

     .. math::
        \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{r}, \tau) =
	 -\Gamma_{A\bar{a}B\bar{b}}(\mathbf{r}, \tau) F_{AB}(\mathbf{r}, \tau) \,,

     where

     .. math::
        F_{AB}(\mathbf{r}, \tau)  =
        \mathcal{F}\big(G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')\big)\,.

     Then we Fourier transform

     .. math::
          \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =
          \mathcal{F}\big(\Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{r}, \tau)\big)\,,

    to get the same result, but with far less computational effort.