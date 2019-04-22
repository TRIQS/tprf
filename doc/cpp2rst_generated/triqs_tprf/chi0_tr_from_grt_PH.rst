..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0_tr_from_grt_PH:

chi0_tr_from_grt_PH
===================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_tr_t chi0_tr_from_grt_PH (triqs_tprf::g_tr_cvt g_tr)

Generalized susceptibility imaginary time bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`


Parameters
----------

 * **g_tr**: Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\tau, \mathbf{r})`.



Return value
------------

Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real-space.

Documentation
-------------


  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
     - G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})