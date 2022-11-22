.. _notation:

Response function notation
==========================

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

   
Generalized susceptibility :math:`\chi`
=======================================

The two particle Green's function can be split into one connected component :math:`G^{(2)}_c` and two dissconnected componens

.. math::
   G^{(2)}(\bar{a}b\bar{c}d) =
   G^{(2)}_c(\bar{a}b\bar{c}d) + G(b\bar{a})G(d\bar{c}) - G(d\bar{a})G(b\bar{c})
   .

Bare generalized susceptibility :math:`\chi_0`
----------------------------------------------
   
The disconnected components can be collected in the bare suceptibilities :math:`\chi^{=}_0` and :math:`\chi^{\times}_0` defined as
   
.. math::
   \chi^{=}_0(\bar{a}b\bar{c}d) \equiv + G(b\bar{a})G(d\bar{c})

.. math::
   \chi^{\times}_0(\bar{a}b\bar{c}d) \equiv - G(d\bar{a})G(b\bar{c})

The non-interacting two-particle response of a system is related to one of these bare susceptibilities depending on what pairs of indices are considered to be in-going and out-going.

Reducible vertex function :math:`F`
-----------------------------------

The remaining connected two-particle Green's function :math:`G^{(2)}_c` is related to the fully reducible vertex function :math:`F` by a dressing of four single-particle propagators, i.e., :math:`G^{(2)}_c \equiv G^4 F`. In detail each external leg of :math:`F(a\bar{a}b\bar{b})` is dressed by a single-particle Green's functions :math:`G`, which gives

.. math::
   G^{(2)}_c(\bar{a}b\bar{c}d)
   \equiv
   G(a\bar{a}) G(b\bar{b}) F(a\bar{b}c\bar{d}) G(c\bar{c}) G(d\bar{d})

Full generalized particle-hole susceptibility :math:`\chi`
----------------------------------------------------------
   
The non-interacting particle-hole excitation response is given by :math:`\chi^{\times}_0` when :math:`\bar{c}d` are incoming and :math:`\bar{a}b` are out-going indices, respectively, and by :math:`\chi^{=}_0` for the in- and out-going pairings :math:`\bar{c}b` and :math:`\bar{a}d`.

The interacting response is given by the generalized susceptibilities :math:`\chi^{\times}` and :math:`\chi^{=}` which are given by the two-particle Green's function after subracting the trivial non-correlated propagation of the in-going and out-going pair

.. math::
   \chi^{=} = G^{(2)} - \chi^{\times}_0 = \chi^{=}_0 + G^{(2)}_c

.. math::
   \chi^{\times} = G^{(2)} - \chi^{=}_0 = \chi^{\times}_0 + G^{(2)}_c

Inserting the equation for the connected two-particle Green's function in terms of the reducible vertex :math:`F` and single particle Green's function gives

.. math::
   \chi^{\times}(\bar{a}b\bar{c}d) =
   \chi^{\times}_0(\bar{a}b\bar{c}d)
   +
   G(a\bar{a})
   G(b\bar{b})
   F(a\bar{b}c\bar{d})
   G(c\bar{c})
   G(d\bar{d})
   \\ =
   \chi^{\times}_0(\bar{a}b\bar{c}d)
   +
   \chi^{\times}_0(\bar{a}b \bar{b}a)
   F(a\bar{b}c\bar{d})
   \chi^{\times}_0(\bar{c}d \bar{d}c)

This response channel is commonly called the Particle-Hole channel (:math:`PH`).

The other possible pairing is given by
   
.. math::
   \chi^{=}(\bar{a}b\bar{c}d) = \chi^{=}_0(\bar{a}b\bar{c}d)
   +
   G(a\bar{a})
   G(d\bar{d})
   F(a\bar{b}c\bar{d})
   G(b\bar{b})
   G(c\bar{c})
   \\ =
   \chi^{=}_0(\bar{a}b\bar{c}d)
   +
   \chi^{=}_0(\bar{a}a\bar{d}d)
   F(a\bar{b}c\bar{d})
   \chi^{=}_0(\bar{b}b \bar{c}c)

and is commonly called the vertical Particle-Hole channel (:math:`\bar{PH}`)
   
Fix me
   
.. math::
   \chi^{\times}(\bar{a}b\bar{c}d) = \chi^{\times}_0(\bar{a}b\bar{c}d)
   + \chi^{\times}_0(\bar{a}b\bar{p}q) F(q\bar{p}r\bar{s}) \chi^{\times}_0(\bar{s}r\bar{c}d)
   
In abbreviated form we can write this in terms of the crossed and direct product :math:`\stackrel{r}{*}` with :math:`r \in \{ =, \times \}`

.. math::
   \chi^{r} = \chi^{r}_0 + \chi^{r}_0 \stackrel{r}{*} F \stackrel{r}{*} \chi^{r}_0

Which we write without explicitly writing down the product channel :math:`r` as

.. math::
   \chi = \chi_0 + \chi_0 F \chi_0
   

Bethe-Salpeter equations (BSE)
==============================
   
The vertex BSEs defines :math:`\chi^r_0` and :math:`\stackrel{r}{*}`

.. math::
   F = \Gamma^r + \Gamma^r \stackrel{r}{*} \chi^r_0 \stackrel{r}{*} F

   

.. math::
   \chi^r = \chi^r_0 + \chi^r_0 \stackrel{r}{*} F \stackrel{r}{*} \chi^r_0


.. math::
   \chi^r = \chi^r_0 + \chi^r_0 \stackrel{r}{*} \Gamma^r \stackrel{r}{*} \chi

.. math::
   \chi_{\bar{a} b \bar{c} d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   \equiv 
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   - G_{b\bar{a}}(\tau_2, \tau_1) G_{d\bar{c}}(\tau_4, \tau_3)

.. math::
   \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4) 
   \equiv
   - G_{d\bar{a}}(\tau_4, \tau_1) G_{b\bar{c}}(\tau_2, \tau_3)

     
Matsubara frequency transforms
------------------------------

Operators and response functions in imaginary time :math:`\tau` can be Fourier transformed to imaginary Matsubara frequencies

.. math::
   \nu_n = \frac{\pi}{\beta}(2n + \vartheta)
   
with :math:`\vartheta = (1-\xi)/2` and :math:`\xi = \pm 1` for bosons/fermions) exploiting :math:`\beta` (anti)periodicity.

The second quantized operators transforms according to

.. math::
   c(i\nu_n) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{i\nu_n \tau} c(\tau)
   \, , \quad
   c^\dagger(i\nu_n) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{-i\nu_n \tau} c^\dagger(\tau)

.. math::
   c(\tau) = \frac{1}{\sqrt{\beta}} \sum_{n=-\infty}^{\infty} e^{-i\nu_n \tau} c(i\nu_n)
   \, , \quad
   c^\dagger(\tau) = \frac{1}{\sqrt{\beta}} \sum_{n=-\infty}^{\infty} e^{i\nu_n \tau} c^\dagger(i\nu_n)

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
   \\ \times
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_1, \tau_2, \tau_3, \tau_4)

   

Particle-hole channel (:math:`PH`)
----------------------------------

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
   \beta \delta_{0, \omega} G_{b\bar{a}}(\nu) G_{d\bar{c}}(\nu')
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
     
Crossed-Particle-particle channel (:math:`PPx`)
-----------------------------------------------

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
   :label: bare_pp_sus_def

.. math::
   \chi^{pp}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   =
   G^{(2), pp}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   - \beta \delta_{\nu+\nu', \omega} G_{b\bar{a}}(\nu) G_{d\bar{c}}(\nu')
