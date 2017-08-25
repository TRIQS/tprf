.. _notation:

Response function notations
===========================

The notation for two-particle response functions and derived quantities such as vertex and the irreducible vertices is not settled and there are many possible notational choices. Here we establish the notation used in TPRF.

The single-particle Green's function :math:`G` is defined as

.. math::
   G_{a\bar{b}}(\tau_1, \tau_2) 
   \equiv 
   - \langle \mathcal{T} c_{a}(\tau_1) c^\dagger_{\bar{b}}(\tau_2) \rangle

and the two-particle Green's function :math:`G^{(2)}` is defined as

.. math::
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   \equiv 
   \langle \mathcal{T} 
   c^\dagger_{\bar{a}}(\tau_1) c_{b}(\tau_2)
   c^\dagger_{\bar{c}}(\tau_3) c_{d}(\tau_4)
   \rangle

The disconnected Green's function :math:`G^{(2),diss}` is given by

.. math::
   G^{(2),diss}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   \equiv
   G_{b\bar{a}}(\tau_2, \tau_1) G_{d\bar{c}}(\tau_4, \tau_3)
   -
   G_{d\bar{a}}(\tau_4, \tau_1) G_{b\bar{c}}(\tau_2, \tau_3)


Abbreviated notation
--------------------

In order to economically express relations invovling two-particle objects we will also employ the short hand notation

.. math::
   G^{(2)}(\bar{a}a\bar{b}b) \equiv
   G^{(2)}_{\bar{a} a \bar{b} b}(\tau_{\bar{a}}, \tau_a, \tau_{\bar{b}}, \tau_b)

and use Einstein summation for repeated labels, e.g.

.. math::
   G(a\bar{a})G(b\bar{a}) \equiv
   \sum_{\bar{a}} \int_0^\beta d \tau_{\bar{a}} \,
   G_{a \bar{a}}(\tau_a, \tau_{\bar{a}}) G_{b \bar{a}}(\tau_b, \tau_{\bar{a}})

Matsubara frequency transforms
------------------------------

Operators and response functions in imaginary time $\tau$ can be Fourier transformed to imagniary Matsubara frequencies :math:`\omega_n = \frac{\pi}{\beta}(2n + \vartheta)` (with :math:`\vartheta = (1-\xi)/2` and :math:`\xi = \pm 1` for bosons/fermions) exploiting :math:`\beta` (anti)periodicity.

The second quantized operators transforms according to

.. math::
   c(i\nu) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{i\nu \tau} c(\tau)
   \, , \quad
   c^\dagger(i\nu) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{-i\nu \tau} c^\dagger(\tau)

The single-particle Green's function :math:`G` transforms as

.. math::
   G_{a\bar{b}}(\nu_1, \nu_2) = \delta_{\nu_1, \nu_2} G_{a\bar{b}}(\nu_1)
   \equiv
   \int_0^\beta d\tau_1 d\tau_2 \,
   \exp \left( i\nu_1 \tau_1 - i \nu_2 \tau_2 \right)
   G_{a\bar{b}}(\tau_1, \tau_2)

The two-particle Green's function :math:`G^{(2)}` transforms according to

.. math::
   G^{(2)}_{\bar{a}b\bar{c}d}(\nu_1, \nu_2, \nu_3, \nu_4)
   =
   \delta_{\nu_1 + \nu_3, \nu_2 + \nu_4}
   G^{(2)}_{\bar{a}b\bar{c}d}(\nu_1, \nu_2, \nu_3, \nu_4)
   \\ \equiv 
   \int_0^\beta d\tau_1 d\tau_2 d\tau_3 d\tau_4
   \exp\left( -i\nu_1 \tau_1 + i \nu_2 \tau_2 - i\nu_3 \tau_3 + i \nu_4 \tau_4 \right)
   \times
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4)

   

Generalized susceptibilities
----------------------------

.. math::
   \chi_{\bar{a} b \bar{c} d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   \equiv 
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   - G_{b\bar{a}}(\tau_2, \tau_1) G_{d\bar{c}}(\tau_4, \tau_3)

.. math::
   \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   \equiv
   - G_{d\bar{a}}(\tau_4, \tau_1) G_{b\bar{c}}(\tau_2, \tau_3)


Particle-hole channel (PH)

.. math::
   \nu_1 = \nu 
   \, , \quad
   \nu_2 = \omega + \nu
   \, , \quad
   \nu_3 = \omega + \nu'
   \, , \quad
   \nu_4 = \nu'

.. math::
   G^{(2),ph}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   =
   G^{(2)}_{\bar{a}b\bar{c}d}(\nu, \omega + \nu, \omega + \nu', \nu')

.. math::
   G^{(2),ph,diss}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   =
   \beta \delta_{\nu+\nu', \omega} G_{b\bar{a}}(\nu) G_{b\bar{c}}(\nu')
   - \beta \delta_{\nu, \nu'} G_{d\bar{a}}(\nu) G_{b\bar{c}}(\omega + \nu)

.. math::
   \chi^{(0),ph}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   =
   - \beta \delta_{\nu, \nu'} G_{d\bar{a}}(\nu) G_{b\bar{c}}(\omega + \nu)

.. math::
   \chi^{ph}_{\bar{a}b\bar{c}d} (\omega, \nu, \nu') 
   =
   G^{(2),ph}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   - \beta \delta_{0, \omega} G_{b\bar{a}}(\nu) G_{d\bar{c}}(\nu')
     
Particle-particle channel (PP)

.. math::
   \nu_1 = \nu
   \, , \quad
   \nu_2 = \omega - \nu'
   \, , \quad
   \nu_3 = \omega - \nu
   \, , \quad
   \nu_4 = \nu'

.. math::
   G^{(2), pp}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') 
   =
   G^{(2)}_{\bar{a}b\bar{c}d}(\nu, \omega - \nu', \omega + \nu, \nu')

.. math::
   G^{(2),pp,diss}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   =
   \beta \delta_{\nu + \nu' , \omega} G_{b\bar{a}}(\nu) G_{d\bar{c}}(\nu')
   - \beta \delta_{\nu, \nu'} G_{d\bar{a}}(\nu) G_{b\bar{c}}(\omega - \nu)

.. math::
   \chi^{(0), pp}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   =
   - \beta \delta_{\nu, \nu'} G_{d\bar{a}}(\nu) G_{b\bar{c}}(\omega - \nu)

.. math::
   \chi^{pp}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   =
   G^{(2), pp}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   - \beta \delta_{\nu+\nu', \omega} G_{b\bar{a}}(\nu) G_{d\bar{c}}(\nu')
