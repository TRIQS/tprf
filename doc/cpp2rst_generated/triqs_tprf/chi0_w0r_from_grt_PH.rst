..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0_w0r_from_grt_PH:

chi0_w0r_from_grt_PH
====================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wr_t chi0_w0r_from_grt_PH (triqs_tprf::g_tr_cvt g_tr)

Generalized susceptibility zero imaginary frequency bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`


Parameters
----------

 * **g_tr**: Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\tau, \mathbf{r})`.



Return value
------------

Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r})` in real-space.

Documentation
-------------


  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r}) =
     - \int_0^\beta d\tau \,
     G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})