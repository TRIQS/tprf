..
   Generated automatically by cpp2rst

.. highlight:: c


.. _gamma_PP_singlet:

gamma_PP_singlet
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t gamma_PP_singlet (triqs_tprf::chi_wk_vt chi_c,
   triqs_tprf::chi_wk_vt chi_s, array_view<std::complex<double>, 4> U_c,
   array_view<std::complex<double>, 4> U_s)

Gamma particle-particle singlet


Parameters
----------

 * **chi_c**: charge susceptibility  :math:`\chi^{(\mathrm{c})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

 * **chi_s**: spin susceptibility  :math:`\chi^{(\mathrm{s})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

 * **U_c**: charge interaction  :math:`U^{(\mathrm{c})}_{a\bar{b}c\bar{d}}`

 * **U_s**: spin interaction  :math:`U^{(\mathrm{s})}_{a\bar{b}c\bar{d}}`



Return value
------------

:math:`\Gamma^{(\mathrm{singlet})}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\omega_n)`

Documentation
-------------


     Computes the particle-particle vertex for singlet pairing in the RPA limit

    .. math::
        \Gamma^{(\mathrm{singlet})}(a\bar{b}c\bar{d}) =
        \frac{3}{2} U^{(\mathrm{s})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{s})}(\bar{B}A\bar{C}D)
        U^{(\mathrm{s})}(D\bar{C}c\bar{d}) \\
        -\frac{1}{2} U^{(\mathrm{c})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{c})}(\bar{B}A\bar{C}D)
        U^{(\mathrm{c})}(D\bar{C}c\bar{d}) \\
       + \frac{1}{2}\big(U^{(\mathrm{s})}(a\bar{b}c\bar{d})+
        U^{(\mathrm{c})}(a\bar{b}c\bar{d})\big)