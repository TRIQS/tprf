..
   Generated automatically by cpp2rst

.. highlight:: c


.. _gw_sigma_tr:

gw_sigma_tr
===========

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_tr_t gw_sigma_tr (triqs_tprf::chi_tr_cvt Wr_tr, triqs_tprf::g_tr_cvt
   g_tr)

GW self energy :math:`\Sigma(\tau, \mathbf{r})` calculator


Parameters
----------

 * **Wr_tr**: retarded screened interaction :math:`W^{(r)}_{abcd}(\tau, \mathbf{r})`

 * **g_tr**: single particle Green's function :math:`G_{ab}(\tau, \mathbf{r})`



Return value
------------

GW self-energy :math:`\Sigma_{ab}(\tau, \mathbf{r})`

Documentation
-------------


    Computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W^{(r)}_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})