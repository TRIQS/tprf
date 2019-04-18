..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_from_gg2_PP:

chi_from_gg2_PP
===============

**Synopsis**:

.. code-block:: c

    triqs_tprf::g2_iw_t chi_from_gg2_PP (triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2)

Generalized susceptibility :math:`\chi^{(0)} = G^{(2)} - GG` in the


Parameters
----------

 * **g**: single particle Green's function :math:`G_{ab}(\nu)`

 * **g2**: two-particle Green's function
   :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`



Return value
------------

chi generalized particle-hole susceptibility
   :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`

Documentation
-------------

   Particle-Particle channel

     Computes

     .. math::
         \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
         - \beta \delta_{\nu + \nu' - \omega} G_{ba}(\nu) \cdot G_{dc}(\nu')