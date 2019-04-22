..
   Generated automatically by cpp2rst

.. highlight:: c


.. _lindhard_chi00_wk:

lindhard_chi00_wk
=================

**Synopsis**:

.. code-block:: c

    triqs_tprf::chi_wk_t lindhard_chi00_wk (triqs_tprf::e_k_cvt e_k, int nw, double beta,
   double mu)

Generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`.


Documentation
-------------


   Analytic calculation of the generalized (non-interacting) Lindhard susceptibility
   in the particle-hole channel. The analytic expression is obtained using residue calculus
   to explicitly evaluate the matsubara sum of the fourier transformed imaginary time
   bubble product of two non-interacting single-particle Green's functions.

   .. math::
      G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n) =
      \left[ i\omega_n \cdot \mathbf{1} - \epsilon(\mathbf{k}) \right]^{-1} .

   The analytic evaluation of the bubble diagram gives

   .. math::
        \chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q}) \equiv
        \mathcal{F} \left\{
	  - G^{(0)}_{d\bar{a}}(\tau, \mathbf{r}) G^{(0)}_{b\bar{c}}(-\tau, -\mathbf{r})
	\right\}
        =
	- \frac{1}{N_k} \sum_{\nu} \sum_{\mathbf{k}}
          G^{(0)}_{d\bar{a}}(\nu, \mathbf{k})
	  G^{(0)}_{b\bar{c}}(\nu + \omega, \mathbf{k} + \mathbf{q})
	\\ =
	- \frac{1}{N_k} \sum_{\nu} \sum_{\mathbf{k}}
	  \left( \sum_{i}
          U^\dagger_{di}(\mathbf{k}) \frac{1}{i\nu - \epsilon_{\mathbf{k}, i}} U_{i\bar{a}}(\mathbf{k})
	  \right)
	  \left( \sum_j
	  U^\dagger_{bj}(\mathbf{k} + \mathbf{q})
	  \frac{1}{i\nu + i\omega - \epsilon_{\mathbf{k} + \mathbf{q}, j}}
	  U_{j\bar{c}}(\mathbf{k} + \mathbf{q})
	  \right)
	\\ =
	\frac{1}{N_k} \sum_{\mathbf{k}} \sum_{ij}
	  \left(
	    [1-\delta(\epsilon_{\mathbf{k},i} - \epsilon_{\mathbf{k}+\mathbf{q}, j})]
	    \frac{ f(\epsilon_{\mathbf{k}, i}) - f(\epsilon_{\mathbf{k}+\mathbf{q}, j}) }
	         {i\omega_n + \epsilon_{\mathbf{k} + \mathbf{q}, j} - \epsilon_{\mathbf{k}, i}}
	    +
  	    \delta(\epsilon_{\mathbf{k},i} - \epsilon_{\mathbf{k}+\mathbf{q}, j})
	    \frac{\beta}{4 \cosh^2 (\beta \epsilon_{\mathbf{k}, i} / 2) }
	  \right)
	  \\ \times
	  U_{i\bar{a}}(\mathbf{k}) U^\dagger_{di}(\mathbf{k})
	  U_{j\bar{c}}(\mathbf{k} + \mathbf{q}) U^\dagger_{bj}(\mathbf{k} + \mathbf{q})

   where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued
   dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

   .. math::
      \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
      = \delta_{ij} \epsilon_{\mathbf{k}, i}

   .. note::
      The analytic formula is sub-optimal in terms of performance for higher temperatures. The evaluation
      scales as :math:`\mathcal{O}(N_k^2)` which is worse than computing the bubble explicitly in imaginary
      time, with scaling :math:`\mathcal{O}(N_k N_\tau \log(N_k N_\tau)` for :math:`N_k \gg N_\tau`.

   .. note::
      Care must be taken when evaluating the fermionic Matsubara frequency sum of the
      product of two simple poles. By extending the sum to an integral over the complex
      plane the standard expression for the Lindhard response is obtained when the
      poles are non-degenerate. The degenerate case produces an additional frequency independent
      contribution (the last term on the last row).