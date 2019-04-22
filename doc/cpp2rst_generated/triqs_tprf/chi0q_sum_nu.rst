..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0q_sum_nu:

chi0q_sum_nu
============

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t chi0q_sum_nu (triqs_tprf::chi_wnk_cvt chi_wnk)

Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})`. (NB! without tail corrections)


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
     \frac{1}{\beta^2} \sum_{\nu=\nu_{min}}^\nu_{max} \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})