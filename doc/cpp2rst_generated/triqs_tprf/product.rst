..
   Generated automatically by cpp2rst

.. highlight:: c


.. _product:

product
=======

**Synopsis**:

.. code-block:: c

    template<Channel_t CH>
    triqs_tprf::g2_iw_t product (triqs_tprf::g2_iw_cvt A, triqs_tprf::g2_iw_cvt B)

Two-particle response-function product :math:`A * B`

.. note::
    Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Template parameters
-------------------

 * **CH**: selects the two-particle channel



Parameters
----------

 * **A**: two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`

 * **B**: two-particle response function :math:`B \equiv A_{abcd}(\omega, \nu, \nu')`



Return value
------------

:math:`(A * B)` in the given channel


Documentation
-------------


 The two-particle response functions :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`
 and :math:`B \equiv B_{abcd}(\omega, \nu, \nu')` are cast to matrix form and their
 product is computed

 .. math::
   (A * B)_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   = \sum_{\bar{\nu}ab}
   A_{\{\nu\alpha\beta\}, \{\bar{\nu}ab\}}(\omega)
   B_{\{\bar{\nu}ab\}, \{\nu'\gamma\delta\}}(\omega)

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the product is returned by value.