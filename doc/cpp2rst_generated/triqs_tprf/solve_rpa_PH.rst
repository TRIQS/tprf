..
   Generated automatically by cpp2rst

.. highlight:: c


.. _solve_rpa_PH:

solve_rpa_PH
============

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t solve_rpa_PH (triqs_tprf::chi_wk_vt chi0,
   array_view<std::complex<double>, 4> U)

Random Phase Approximation (RPA) in the particle-hole channel


Parameters
----------

 * **chi0**: bare particle-hole bubble :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

 * **U**: RPA static vertex as obtained from triqs_tprf.rpa_tensor.get_rpa_tensor :math:`U_{a\bar{b}c\bar{d}}`



Return value
------------

RPA suceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

Documentation
-------------


     Computes the equation

     .. math::
         \chi(\bar{a}b\bar{c}d) = \big(
         \mathbb{1}
         - \chi^{(0)}(\bar{a}b\bar{B}A) U(A\bar{B}D\bar{C})
         \big)^{-1} \chi^{(0)}(\bar{C}D\bar{c}d)\,.