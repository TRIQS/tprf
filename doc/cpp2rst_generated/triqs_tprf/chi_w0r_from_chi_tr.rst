..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_w0r_from_chi_tr:

chi_w0r_from_chi_tr
===================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wr_t chi_w0r_from_chi_tr (triqs_tprf::chi_tr_vt chi_tr)

Static susceptibility calculation :math:`\chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`


Parameters
----------

 * **chi_tr**: Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real space



Return value
------------

Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real-space

Documentation
-------------


  Explicit calculation of the static, zero frequency response, by 2nd order trapetzoidal
  integration in imaginary time, i.e.

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) =
         \int_0^\beta d\tau \, \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})