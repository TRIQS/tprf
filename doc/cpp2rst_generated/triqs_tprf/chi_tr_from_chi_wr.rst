..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_tr_from_chi_wr:

chi_tr_from_chi_wr
==================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_tr_t chi_tr_from_chi_wr (triqs_tprf::chi_wr_cvt chi_wr, int ntau = -1)

Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`


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
         \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
         \mathcal{F}_{\omega \rightarrow \tau} \left\{
	 \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \right\}