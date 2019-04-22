..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_wr_from_chi_wk:

chi_wr_from_chi_wk
==================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wr_t chi_wr_from_chi_wk (triqs_tprf::chi_wk_cvt chi_wk)

Parallell Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`


Parameters
----------

 * **chi_wr**: Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`
                in imaginary time and momentum space.



Return value
------------

Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`
          in Matsubara frequency and real space.

Documentation
-------------


  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \mathcal{F}_{\mathbf{k} \rightarrow \mathbf{r}} \left\{
         \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})
         \right\}