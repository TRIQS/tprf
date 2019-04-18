..
   Generated automatically by cpp2rst

.. highlight:: c


.. _gw_self_energy:

gw_self_energy
==============

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_wk_t gw_self_energy (triqs_tprf::chi_wk_vt W_wk, triqs_tprf::g_wk_vt
   g_wk)

GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator


Parameters
----------

 * **W_wk**: screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`

 * **g_wk**: single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`



Return value
------------

GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`

Documentation
-------------


    Fourier transforms the screened interaction and the single-particle
    Green's function to imagiary time and real space.

    .. math::
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

    .. math::
        W_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W_{abcd}(i\omega_n, \mathbf{k}) \right\}

    computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    and transforms back to frequency and momentum

    .. math::
        \Sigma_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma_{ab}(\tau, \mathbf{r}) \right\}