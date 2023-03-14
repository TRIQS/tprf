/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include "../types.hpp"

namespace triqs_tprf {

  /** Generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`.
   
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
             [1 - \delta_{0, \omega_n} \delta_{\epsilon_{\mathbf{k},i},\epsilon_{\mathbf{k}+\mathbf{q}, j}})]
             \frac{ f(\epsilon_{\mathbf{k}, i}) - f(\epsilon_{\mathbf{k}+\mathbf{q}, j}) }
                  {i\omega_n + \epsilon_{\mathbf{k} + \mathbf{q}, j} - \epsilon_{\mathbf{k}, i}}
             +
             \delta_{0, \omega_n} \delta_{\epsilon_{\mathbf{k},i},\epsilon_{\mathbf{k}+\mathbf{q}, j}}
             \frac{\beta}{4 \cosh^2 (\beta \epsilon_{\mathbf{k}, i} / 2) }
           \right)
           \\ \times
           U_{\bar{a}i}(\mathbf{k}) U^\dagger_{id}(\mathbf{k}) 
           U_{\bar{c}j}(\mathbf{k} + \mathbf{q}) U^\dagger_{jb}(\mathbf{k} + \mathbf{q})

    where the $U(\mathbf{k})$ matrices are the diagonalizing unitary transform of the matrix valued 
    dispersion relation $\epsilon_{\bar{a}b}(\mathbf{k})$, i.e.

    @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
    @param mesh bosonic Matsubara frequency mesh
    @param mu chemical potential :math:`\mu`
    @return generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`

    .. math::
       \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}

    .. note::
       The analytic formula is sub-optimal in terms of performance for higher temperatures. The evaluation
       scales as $\mathcal{O}(N_k^2)$ which is worse than computing the bubble explicitly in imaginary 
       time, with scaling $\mathcal{O}(N_k N_\tau \log(N_k N_\tau)$ for $N_k \gg N_\tau$.

    .. note::
       Care must be taken when evaluating the fermionic Matsubara frequency sum of the
       product of two simple poles. By extending the sum to an integral over the complex 
       plane the standard expression for the Lindhard response is obtained when the 
       poles are non-degenerate. The degenerate case produces an additional frequency independent
       contribution (the last term on the last row).
*/
  chi_wk_t lindhard_chi00(e_k_cvt e_k, gf_mesh<imfreq> mesh, double mu);

  /** Generalized Lindhard susceptibility in the particle-hole channel and for real frequencies :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q})`.
    
    Analytic calculation of the generalized (non-interacting) Lindhard susceptibility 
    in the particle-hole channel in real frequencies. The analytic expression is obtained using 
    residue calculus to explicitly evaluate the matsubara sum of the fourier transformed imaginary
    time bubble product of two non-interacting single-particle Green's functions.

    .. math::
       G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n) =
       \left[ i\omega_n \cdot \mathbf{1} - \epsilon(\mathbf{k}) \right]^{-1} .

    The analytic continuation of the resulting expression to the real frequency axis gives
 
    .. math::
        \chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q}) = 
        \frac{1}{N_k} \sum_{\mathbf{k}} \sum_{ij} 
        \frac{ f(\epsilon_{\mathbf{k}, i}) - f(\epsilon_{\mathbf{k}+\mathbf{q}, j}) }
            {\omega + i\delta + \epsilon_{\mathbf{k} + \mathbf{q}, j} - \epsilon_{\mathbf{k}, i}}
        \\ \times
        U_{\bar{a}i}(\mathbf{k}) U^\dagger_{id}(\mathbf{k}) 
        U_{\bar{c}j}(\mathbf{k} + \mathbf{q}) U^\dagger_{jb}(\mathbf{k} + \mathbf{q})
 
    where the $U(\mathbf{k})$ matrices are the diagonalizing unitary transform of the matrix valued 
    dispersion relation $\epsilon_{\bar{a}b}(\mathbf{k})$, i.e.
 
    .. math::
       \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}

    @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
    @param mesh real frequency mesh 
    @param beta inverse temperature
    @param mu chemical potential :math:`\mu`
    @param delta broadening :math:`\delta`
    @return real frequency generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q})`
*/
  chi_fk_t lindhard_chi00(e_k_cvt e_k, gf_mesh<refreq> mesh, double beta, double mu, double delta);

} // namespace triqs_tprf
