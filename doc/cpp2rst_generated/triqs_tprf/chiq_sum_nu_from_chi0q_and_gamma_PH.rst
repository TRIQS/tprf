..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chiq_sum_nu_from_chi0q_and_gamma_PH:

chiq_sum_nu_from_chi0q_and_gamma_PH
===================================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_kw_t chiq_sum_nu_from_chi0q_and_gamma_PH (triqs_tprf::chi_wnk_cvt
   chi0_wnk, triqs_tprf::chi_wnn_cvt gamma_ph_wnn)

Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.


Parameters
----------

 * **chi0_wnk**: Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

 * **gamma_ph_wnn**: Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.



Return value
------------

Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

Documentation
-------------


  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \chi^{(0)} \left[ 1 - \Gamma^{(PH)} \chi^{(0)} \right]^{-1}