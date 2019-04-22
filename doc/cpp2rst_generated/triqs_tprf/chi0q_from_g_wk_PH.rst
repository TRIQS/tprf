..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi0q_from_g_wk_PH:

chi0q_from_g_wk_PH
==================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wnk_t chi0q_from_g_wk_PH (int nw, int nn, triqs_tprf::g_wk_cvt g_wk)

Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` with convolution in k-space.


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
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
     - \frac{\beta}{N_k} \sum_\mathbf{k}
     G_{d\bar{a}}(\nu, \mathbf{k}) \cdot G_{b\bar{c}}(\nu + \omega, \mathbf{k} - \mathbf{q})