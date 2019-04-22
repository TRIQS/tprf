..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_wr_from_chi_tr:

chi_wr_from_chi_tr
==================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wr_t chi_wr_from_chi_tr (triqs_tprf::chi_tr_cvt chi_tr, int nw)

Parallell Fourier transform from  :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`


Parameters
----------

 * **chi_tr**: Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`
                in imaginary time and real space.



Return value
------------

Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`
          in Matsubara frequency and real-space.

Documentation
-------------


  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \mathcal{F}_{\tau \rightarrow \omega} \left\{
         \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})
         \right\}