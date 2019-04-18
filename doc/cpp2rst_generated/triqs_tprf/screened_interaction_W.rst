..
   Generated automatically by cpp2rst

.. highlight:: c


.. _screened_interaction_W:

screened_interaction_W
======================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t screened_interaction_W (triqs_tprf::chi_wk_vt PI_wk,
   triqs_tprf::chi_k_vt V_k)

Screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator


Parameters
----------

 * **PI_wk**: polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`

 * **V_k**: static interaction :math:`V_{abcd}(\mathbf{k})`



Return value
------------

screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`

Documentation
-------------

    for static momentum dependent interaction :math:`V(\mathbf{k})`

    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) =
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          W_{hgcd}(i\omega_n, \mathbf{k})