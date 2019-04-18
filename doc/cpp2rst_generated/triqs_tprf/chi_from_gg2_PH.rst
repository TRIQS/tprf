..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_from_gg2_PH:

chi_from_gg2_PH
===============

**Synopsis**:

.. code-block:: c

    triqs_tprf::g2_iw_t chi_from_gg2_PH (triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2)

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

   Particle-Hole channel

     Computes

     .. math::
         \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
         - \beta \delta_{\omega} G_{ba}(\nu) \cdot G_{dc}(\nu')