..
   Generated automatically by cpp2rst

.. highlight:: c


.. _inverse:

inverse
=======

**Synopsis**:

.. code-block:: c

    template<Channel_t CH>
    triqs_tprf::g2_iw_t inverse (triqs_tprf::g2_iw_cvt g)

Two-particle response-function inversion $[g]^{-1}$.

.. note::
    Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Template parameters
-------------------

 * **CH**: selects the two-particle channel



Parameters
----------

 * **g**: two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`



Return value
------------

:math:`[g]^{-1}` in the given channel


Documentation
-------------


 The two-particle response function :math:`g_{abcd}(\omega, \nu, \nu')`
 is cast to matrix form and inverted

 .. math::
   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the inverse is returned by value.