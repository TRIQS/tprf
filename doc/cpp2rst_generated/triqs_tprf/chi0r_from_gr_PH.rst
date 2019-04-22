..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0r_from_gr_PH:

chi0r_from_gr_PH
================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wnr_t chi0r_from_gr_PH (int nw, int nn, triqs_tprf::g_wr_cvt g_nr)

Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})`.


Parameters
----------

 * **nw**: Number of bosonic Matsubara freqiencies.

 * **nn**: Number of fermionic Matsubara freqiencies.

 * **g_tr**: Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.



Return value
------------

Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.

Documentation
-------------


  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     - \beta G_{d\bar{a}}(\nu, \mathbf{r}) \cdot G_{b\bar{c}}(\nu + \omega, -\mathbf{r})