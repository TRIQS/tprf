..
   Generated automatically by cpp2rst

.. highlight:: c


.. _dynamical_screened_interaction_W_wk_from_generalized_susceptibility:

dynamical_screened_interaction_W_wk_from_generalized_susceptibility
===================================================================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t
   dynamical_screened_interaction_W_wk_from_generalized_susceptibility
   (triqs_tprf::chi_wk_cvt chi_wk, triqs_tprf::chi_k_cvt V_k)

Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator


Parameters
----------

 * **chi_wk**: polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`

 * **V_k**: static interaction :math:`V_{abcd}(\mathbf{k})`



Return value
------------

dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`

Documentation
-------------

    for static momentum-dependent interactions :math:`V(\mathbf{k})` and
    known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          V_{hgcd}(\mathbf{k})

    Instead of returning :math:`W^{(full)}` we return the dynamical/retarded part
    :math:`W^{(r)}` (with zero high-frequency offset)

    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) =
            W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) - V_{abcd}(\mathbf{k})