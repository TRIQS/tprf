..
   Generated automatically by cpp2rst

.. highlight:: c


.. _fourier_wr_to_tr:

fourier_wr_to_tr
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_tr_t fourier_wr_to_tr (triqs_tprf::g_wr_cvt g_wr, int nt = -1)

Fast fourier transform of real-space Green's function from Matsubara frequency to imaginary time


Parameters
----------

 * **g_wr**: real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`



Return value
------------

real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`

Documentation
-------------


    Computes: :math:`G_{a\bar{b}}(\tau, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`