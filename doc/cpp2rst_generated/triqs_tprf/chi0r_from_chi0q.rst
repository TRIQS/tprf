..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0r_from_chi0q:

chi0r_from_chi0q
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wnr_t chi0r_from_chi0q (triqs_tprf::chi_wnk_cvt chi_wnk)

Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum-space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real-space.


Parameters
----------

 * **chi_wnk**: Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.



Return value
------------

Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.

Documentation
-------------


  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     \mathcal{F}_{\mathbf{q} \rightarrow \mathbf{r}} \left\{
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})
     \right\}