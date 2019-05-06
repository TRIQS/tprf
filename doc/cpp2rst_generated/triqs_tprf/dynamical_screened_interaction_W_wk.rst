..
   Generated automatically by cpp2rst

.. highlight:: c


.. _dynamical_screened_interaction_W_wk:

dynamical_screened_interaction_W_wk
===================================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t dynamical_screened_interaction_W_wk (triqs_tprf::chi_wk_cvt
   PI_wk, triqs_tprf::chi_k_cvt V_k)

Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator


Parameters
----------

 * **PI_wk**: polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`

 * **V_k**: static interaction :math:`V_{abcd}(\mathbf{k})`



Return value
------------

dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`

Documentation
-------------

    for static momentum-dependent interactions :math:`V(\mathbf{k})`.

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

    Instead of returning :math:`W^{(full)}` we return the dynamical/retarded part
    :math:`W^{(r)}` (with zero high-frequency offset)

    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) =
            W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) - V_{abcd}(\mathbf{k})