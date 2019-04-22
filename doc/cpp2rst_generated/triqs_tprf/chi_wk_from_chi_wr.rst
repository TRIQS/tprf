..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_wk_from_chi_wr:

chi_wk_from_chi_wr
==================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t chi_wk_from_chi_wr (triqs_tprf::chi_wr_cvt chi_wr)

Parallell Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`


Parameters
----------

 * **chi_wr**: Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`
                in Matsubara frequency and real space.



Return value
------------

Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`
          in Matsubara frequency and momentum space.

Documentation
-------------


  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
         \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{k}} \left\{
         \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})
         \right\}