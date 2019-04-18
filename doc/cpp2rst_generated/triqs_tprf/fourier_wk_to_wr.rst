..
   Generated automatically by cpp2rst

.. highlight:: c


.. _fourier_wk_to_wr:

fourier_wk_to_wr
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_wr_t fourier_wk_to_wr (triqs_tprf::g_wk_cvt g_wk)

Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space


Parameters
----------

 * **g_wk**: k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`



Return value
------------

real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`

Documentation
-------------


    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F}^{-1} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \right\}`