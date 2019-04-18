..
   Generated automatically by cpp2rst

.. highlight:: c


.. _eliashberg_product:

eliashberg_product
==================

**Synopsis**:

.. code-block:: c

    triqs_tprf::gk_iw_t eliashberg_product (triqs_tprf::chi_wk_vt Gamma_pp,
   triqs_tprf::gk_iw_vt g_wk, triqs_tprf::gk_iw_vt delta_wk)

Linearized Eliashberg product


Parameters
----------

 * **chi_pp**: particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\nu_n)`

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
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')