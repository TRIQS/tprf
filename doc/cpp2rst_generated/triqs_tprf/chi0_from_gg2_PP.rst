..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0_from_gg2_PP:

chi0_from_gg2_PP
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g2_iw_t chi0_from_gg2_PP (triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt
   g2)

Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Particle


Parameters
----------

 * **g**: single particle Green's function :math:`G_{ab}(\nu)`

 * **g2**: two-particle Green's function with the mesh to use for
   :math:`\chi^{(0)}`



Return value
------------

chi0 particle-particle bubble
   :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`

Documentation
-------------

   channel

     Computes

     .. math::
         \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         - \beta \delta_{\nu, \nu'} G_{da}(\nu) \cdot G_{bc}(\omega - \nu)