..
   Generated automatically by cpp2rst

.. highlight:: c


.. _fourier_tr_to_wr:

fourier_tr_to_wr
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_wr_t fourier_tr_to_wr (triqs_tprf::g_tr_cvt g_tr, int nw = -1)

Fast fourier transform of real-space Green's function from imaginary time to Matsubara frequency


Parameters
----------

 * **g_tr**: real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`



Return value
------------

real-space Matsubara frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`

Documentation
-------------


    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(\tau, \mathbf{r}) \right\}`