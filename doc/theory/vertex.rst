.. _vertex:

Vertex functions
================

.. note:: This notation follows closely [Ayral, Parcollet, PRB 94, 075159 (2016)] with the exception that :math:`PH` and :math:`\bar{PH}` are interchanged.

Fully reducible vertex :math:`F`
--------------------------------

The fully recucible vertex function :math:`F(a\bar{a}b\bar{b})` is defined as the amputation of the connected two-particle Green's function :math:`G^{(2)}_c` by four single-particle Green's functions :math:`G`, one for each leg.

.. math::
   G(a\bar{a}) G(b\bar{b}) F(a\bar{b}c\bar{d}) G(c\bar{c}) G(d\bar{d})
   \equiv G^{(2)}_c(\bar{a}b\bar{c}d)

Bethe-Salpeter equations for the fully reducible vertext :math:`F`
==================================================================
   
The Bethe-Salpeter equations for the fully reducible vertex :math:`F` is defined for any given channel :math:`r` as

.. math::
   F = \Gamma^{r} + \Gamma^{r} GG F

The possible pairings of indices in the right-hand side prodict produces three non-equivalent equations labeled :math:`r \in \{ PH, \bar{PH}, PP\}` standing for, the particle-hole (:math:`PH`), vertical-particle-hole (:math:`\bar{PH}`), and particle-particle (:math:`PP`) channel, respectively.

Each equation and index pairing is associated with one "channel" within which the :math:`r`-channel-irreducible vertex :math:`\Gamma^r`.

.. math::
   F(a\bar{b}c\bar{d}) = \Gamma^{PH}(a\bar{b}c\bar{d}) -
   \Gamma^{PH}(a\bar{b}u\bar{v}) G(u\bar{u}) G(v\bar{v}) F(v\bar{u}c\bar{d})

.. math::
   F(a\bar{b}c\bar{d}) = \Gamma^{\bar{PH}}(a\bar{b}c\bar{d}) +
   \Gamma^{\bar{PH}}(a\bar{u}v\bar{d}) G(u\bar{u}) G(v\bar{v}) F(u\bar{b}c\bar{v})
   
.. math::
   F(a\bar{b}c\bar{d}) = \Gamma^{PP}(a\bar{b}c\bar{d}) + \frac{1}{2}
   \Gamma^{PP}(a\bar{u}c\bar{v}) G(u\bar{u}) G(v\bar{v}) F(v\bar{b}u\bar{d})

Collecting the two single-particle Green's functions into a channel dependent "bare" susceptibility :math:`\chi^r_0` the Behte-Salpeter equations can be expressed as

.. math::
   F = \Gamma^r + \Gamma^r \stackrel{r}{*} \chi^r_0 \stackrel{r}{*}  F 

where the channel non-interacting vertex functions :math:`\chi_0^r` are defined as

.. math::
   \chi_0^{PH}(\bar{a}b\bar{c}d) = G(b\bar{c}) G(d\bar{a})
   \, , \quad
   \big(
   \chi_0^{PH}(\bar{v}u\bar{u}v) = G(u\bar{u}) G(v\bar{v})
   \big)

.. math::
   \chi_0^{\bar{PH}}(\bar{a}b\bar{c}d) = G(b\bar{a}) G(d\bar{c})
   \, , \quad
   \big(
   \chi_0^{\bar{PH}}(\bar{u}u\bar{v}v) = G(u\bar{u}) G(v\bar{v})
   \big)

.. math::
   \chi_0^{PP}(\bar{a}b\bar{c}d) = G(d\bar{a}) G(b\bar{c})
   \, , \quad
   \big(
   \chi_0^{PP}(\bar{u}v\bar{v}u) = G(u\bar{u}) G(v\bar{v})
   \big)

yielding the Bethe-Salpeter equations expressed fully in two-particle quantities

.. math::
   F(a\bar{b}c\bar{d}) = \Gamma^{PH}(a\bar{b}c\bar{d}) -
   \Gamma^{PH}(a\bar{b}q\bar{p}) \, \chi^{PH}_0(\bar{p}q\bar{r}s) \, F(s\bar{r}c\bar{d})
   :label: BSE_PH

.. math::
   F(a\bar{b}c\bar{d}) = \Gamma^{\bar{PH}}(a\bar{b}c\bar{d}) +
   \Gamma^{\bar{PH}}(a\bar{p}s\bar{d}) \, \chi^{\bar{PH}}_0(\bar{p}q\bar{r}s) \, F(q\bar{b}c\bar{r})
   :label: BSE_PHbar
   
.. math::
   F(a\bar{b}c\bar{d}) = \Gamma^{PP}(a\bar{b}c\bar{d}) + \frac{1}{2}
   \Gamma^{PP}(a\bar{p}c\bar{r}) \, \chi^{PP}_0(\bar{p}q\bar{r}s) \, F(q\bar{b}s\bar{d})
   :label: BSE_PP

Note that the notion of the product :math:`\stackrel{r}{*}` is dependent on the channel :math:`r` and consists of different choices of contracting two indices in Eqs. :eq:`BSE_PH`, :eq:`BSE_PHbar`, and :eq:`BSE_PP`.


Matsubara frequency parametrization
===================================

In a carefully choosen parametrization of Matsubara frequencies, the two-time integrals appearing in the products of the Bethe-Salpeter equations :eq:`BSE_PH`, :eq:`BSE_PHbar`, and :eq:`BSE_PP` can be reduced to a single sum over one Matsubara frequency. This is achieved by using a channel dependent three frequency reparametrization that directly imposes total frequency conservation, the forms are

.. math::
   \begin{array}{ll}
   PH: \left\{
   \begin{array}{rl}
   \nu_1 &=& \nu \\
   \nu_2 &=& \nu + \omega \\
   \nu_3 &=& \nu' + \omega \\
   \nu_4 &=& \nu'  
   \end{array}
   \right.
   \, , & \quad
   \bar{PH}: \left\{
   \begin{array}{rcl}
   \nu_1 &=& \nu \\
   \nu_2 &=& \nu'\\
   \nu_3 &=& \nu' + \omega\\
   \nu_4 &=& \nu + \omega 
   \end{array}\right.
   \, , \quad
   \\ \\
   PP: \left\{
   \begin{array}{rcl}
   \nu_1 &=& \nu \\
   \nu_2 &=& \nu' \\
   \nu_3 &=& \omega - \nu \\
   \nu_4 &=& \omega - \nu'
   \end{array}\right.
   \, , & \quad   
   PPx: \left\{
   \begin{array}{rcl}
   \nu_1 &=& \nu \\
   \nu_2 &=& \omega - \nu' \\
   \nu_3 &=& \omega - \nu \\
   \nu_4 &=& \nu'
   \end{array}\right.
   \end{array}
   :label: channel_freq

for the (horizontal) Particle-Hole (:math:`PH`) channel, the (vertical) Particle-Hole (:math:`\bar{PH}`) channel, the Particle-Particle (:math:`PP`) channel, and the Crossed-Particle-Particle (:math:`PPx`) channel, respectively.

.. The (horizontal) Particle-Hole (:math:`PH`) channel:

.. 
   \begin{align}
   \nu_1 &= \nu \, , \\
   \nu_2 &= \nu + \omega \, , \\
   \nu_3 &= \nu' + \omega \, , \\
   \nu_4 &= \nu' \, .
   \end{align}
   :label: PH_freq
   
.. The (vertical) Particle-Hole (:math:`\bar{PH}`) channel:

.. 
   \begin{align}
   \nu_1 &= \nu \, , \\
   \nu_2 &= \nu' \, , \\
   \nu_3 &= \nu' + \omega \, , \\
   \nu_4 &= \nu + \omega \, .
   \end{align}
   :label: PHbar_freq
   
.. The Particle-Particle (:math:`PP`) channel:

.. 
   \begin{align}
   \nu_1 &= \nu \, , \\
   \nu_2 &= \nu' \, , \\
   \nu_3 &= \omega - \nu \, , \\
   \nu_4 &= \omega - \nu' \, .
   \end{align}
   :label: PP_freq

.. note::
   The current definition of the :math:`PP` channel in `cthyb` is actually the :math:`PPx` channel. FIXME!
	   
.. The Crossed-Particle-Particle (:math:`PPx`) channel:

.. 
   \begin{align}
   \nu_1 &= \nu \, , \\
   \nu_2 &= \omega - \nu' \, , \\
   \nu_3 &= \omega - \nu \, , \\
   \nu_4 &= \nu' \, .
   \end{align}
   :label: PPx_freq

In terms of imaginary time the channel dependent three frequency representation maps to the follwing pairing of the four imaginary times :math:`\tau_a`, :math:`\tau_\bar{b}`, :math:`\tau_c`, :math:`tau_{\bar{d}}` of a response function :math:`\chi_{a\bar{b}c\bar{d}}(\tau_a, \tau_{\bar{b}}, \tau_c, \tau_{\bar{d}})`

.. math::
   PH : \,
   + i\omega (\tau_{\bar{b}} - \tau_c)
   + i\nu    (-\tau_{a} + \tau_{\bar{b}})
   + i\nu'   (-\tau_{c} + \tau_{\bar{d}})

.. math::
   \bar{PH} : \,
   + i\omega (-\tau_{c} + \tau_{\bar{d}})
   + i\nu    (-\tau_{a} + \tau_{\bar{d}})
   + i\nu'   (\tau_{\bar{b}} - \tau_{c})

.. math::
   PP : \,
   + i\omega (-\tau_{c} + \tau_{\bar{d}})
   + i\nu    (-\tau_{a} + \tau_{c})
   + i\nu'   (\tau_{\bar{b}} - \tau_{\bar{d}})

.. math::
   PPx : \,
   + i\omega (\tau_{\bar{b}} - \tau_{c})
   + i\nu    (-\tau_{a} + \tau_{c})
   + i\nu'   (-\tau_{\bar{b}} + \tau_{\bar{d}})

In a general product :math:`P = \Gamma \stackrel{r}{*} \chi_0` the total frequency conservation of the components of the product :math:`\Gamma` and :math:`\chi_0` gives two constraints that when combined gives the total frequency conservation of the product :math:`P` and a reduction of the frequency summation of the product from two frequencies to one. This is achieved by using the above global reparametrizations of the four fermionic Matsubara frequencies :math:`\nu_1 ,\, \nu_2 ,\, \nu_3 ,\, \nu_4` of every response function :math:`Q(\nu_1\nu_2\nu_3\nu_4)` for the particular channel :math:`r \in \{PH, \bar{PH}, PP\}` in question.

In order to map the products to matrix products in index and frequency space the following index ordering has to be done

.. math::
   PH: \,
   \chi^{PH}_{a\bar{b}c\bar{d}}(\omega, \nu, nu') =
   \chi^{PH}_{\{\bar{\nu}, \bar{a}b \},\{ \nu, d\bar{c} \}}(\omega)

.. math::
   \bar{PH}: \,
   \chi^{\bar{PH}}_{a\bar{b}c\bar{d}}(\omega, \nu, nu') =
   \chi^{\bar{PH}}_{\{\bar{\nu}, \bar{a}d \}, \{\nu', b\bar{c} \} }(\omega)

.. math::
   PP: \,
   \chi^{PP}_{a\bar{b}c\bar{d}}(\omega, \nu, nu') =
   \chi^{PP}_{\{ \bar{\nu}, \bar{a}\bar{c} \}, \{ \nu', bd\}}
   (\omega)

.. math::
   PPx: \,
   \chi^{PPx}_{a\bar{b}c\bar{d}}(\omega, \nu, nu') =
   \chi^{PPx}_{\{ \bar{\nu}, \bar{c}\bar{a} \}, \{ \nu', bd\}}
   (\omega)

The resulting productformulas reads (see separate derivation chapter),

.. math::
   P^{PH}_{a\bar{b}\bar{c}d}(\omega, \nu,\nu') =
   \frac{1}{\beta^2} \sum_{\bar{\nu} u\bar{v}}
   \Gamma^{PH}_{ \{ \nu, a\bar{b} \},\{ \bar{\nu}, \bar{v}u \}}(\omega)
   \,
   \chi^{PH}_{0, \{\bar{\nu}, \bar{v}u \},\{ \nu, d\bar{c} \}}(\omega)

.. math::
   P^{\bar{PH}}_{ab\bar{c}\bar{d}}(\omega, \nu, \nu')
   =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u} v}
   \Gamma^{\bar{PH}}_{\{ \nu, a\bar{d} \}, \{ \bar{\nu}, \bar{u}v \}}(\omega)
   \,
   \chi^{\bar{PH}}_{0, \{\bar{\nu}, \bar{u}v \}, \{\nu', b\bar{c} \} }(\omega)

.. math::
   P^{PP}_{abcd}(\omega, \nu, \nu')
   =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{ \{ \nu , ac \}, \{\bar{\nu}, \bar{u}\bar{v} \} }
   (\omega)
   \,
   \chi^{PP}_{0, \{ \bar{\nu}, \bar{u}\bar{v} \}, \{ \nu', bd\}}
   (\omega)

.. math::
   P^{PPx}_{abcd}(\omega, \nu, \nu')
   =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PPx}_{ \{ \nu , ca \}, \{\bar{\nu}, \bar{u}\bar{v} \} }
   (\omega)
   \,
   \chi^{PPx}_{0, \{ \bar{\nu}, \bar{u}\bar{v} \}, \{ \nu', bd\}}
   (\omega)

