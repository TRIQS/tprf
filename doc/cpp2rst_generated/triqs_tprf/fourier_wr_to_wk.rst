..
   Generated automatically by cpp2rst

.. highlight:: c


.. _fourier_wr_to_wk:

fourier_wr_to_wk
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_wk_t fourier_wr_to_wk (triqs_tprf::g_wr_cvt g_wr)

Fast fourier transform of imaginary frequency Green's function from real-space to k-space


Parameters
----------

 * **g_wr**: real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`



Return value
------------

k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

Documentation
-------------


    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`