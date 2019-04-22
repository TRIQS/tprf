..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0q_from_chi0r:

chi0q_from_chi0r
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wnk_t chi0q_from_chi0r (triqs_tprf::chi_wnr_cvt chi_wnr)

Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum space.


Parameters
----------

 * **chi_wnr**: Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.



Return value
------------

Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Documentation
-------------


  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
     \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{q}} \left\{
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})
     \right\}