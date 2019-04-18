..
   Generated automatically by cpp2rst

.. highlight:: c


.. _identity:

identity
========

**Synopsis**:

.. code-block:: c

    template<Channel_t CH>
    triqs_tprf::g2_iw_t identity (triqs_tprf::g2_iw_cvt g)

Two-particle response-function identity operator :math:`\mathbf{1}`

.. note::
    Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Template parameters
-------------------

 * **CH**: selects the two-particle channel



Parameters
----------

 * **A**: two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')` determinig the shape and size of the unity operator



Return value
------------

the unity operator :math:`\mathbf{1}`, in the given channel


Documentation
-------------


 Constructs the unity-operator in the given channel

 .. math::
   \mathbf{1}_{abcd}(\omega,\nu,\nu') =
   \mathbf{1}_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   \equiv
   \delta_{\nu\nu'} \delta_{\alpha\gamma} \delta_{\beta\delta}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the result is returned by value.