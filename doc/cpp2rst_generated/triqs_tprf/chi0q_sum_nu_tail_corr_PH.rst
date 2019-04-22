..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0q_sum_nu_tail_corr_PH:

chi0q_sum_nu_tail_corr_PH
=========================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t chi0q_sum_nu_tail_corr_PH (triqs_tprf::chi_wnk_cvt chi_wnk)

Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})` using higher order tail corrections when summing to infinity.


Parameters
----------

 * **chi_wnk**: Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.



Return value
------------

Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic Matsubara frequency and momentum space.

Documentation
-------------


  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \frac{1}{\beta^2} \sum_{\nu=-\infty}^\infty \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})