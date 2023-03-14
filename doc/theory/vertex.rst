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

Bethe-Salpeter equations for the fully reducible vertex :math:`F`
-----------------------------------------------------------------
  
The Bethe-Salpeter equations for the fully reducible vertex :math:`F` is defined for any given channel :math:`r` as

.. math::
   F = \Gamma^{r} + \Gamma^{r} \chi_0^r F

The possible pairings of indices in the right-hand side product produces three non-equivalent equations labeled :math:`r \in \{ PH, \bar{PH}, PP\}` standing for, the particle-hole (:math:`PH`), vertical-particle-hole (:math:`\bar{PH}`), and particle-particle (:math:`PP`) channel, respectively.

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

Collecting the two single-particle Green's functions into a channel dependent "bare" susceptibility :math:`\chi^r_0` the Bethe-Salpeter equations can be expressed as

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
-----------------------------------

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
   :label: matsubara_vertex

for the (horizontal) Particle-Hole (:math:`PH`) channel, the (vertical) Particle-Hole (:math:`\bar{PH}`) channel, the Particle-Particle (:math:`PP`) channel, and the Crossed-Particle-Particle (:math:`PPx`) channel, respectively.


In terms of imaginary time the channel dependent three frequency representation maps to the follwing pairing of the four imaginary times :math:`\tau_a`, :math:`\tau_\bar{b}`, :math:`\tau_c`, :math:`\tau_{\bar{d}}` of a response function :math:`\chi_{a\bar{b}c\bar{d}}(\tau_a, \tau_{\bar{b}}, \tau_c, \tau_{\bar{d}})`

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

In a general product :math:`P = \Gamma \stackrel{r}{*} \chi_0`, the total frequency conservation of the components of the product :math:`\Gamma` and :math:`\chi_0` gives two constraints that when combined gives the total frequency conservation of the product :math:`P` and a reduction of the frequency summation of the product from two frequencies to one. This is achieved by using the above global reparametrizations of the four fermionic Matsubara frequencies :math:`\nu_1 ,\, \nu_2 ,\, \nu_3 ,\, \nu_4` of every response function :math:`Q(\nu_1\nu_2\nu_3\nu_4)` for the particular channel :math:`r \in \{PH, \bar{PH}, PP\}` in question.

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

The resulting product formulas reads (see separate derivation chapter),

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


.. _derivation_index_pairing:

Derivation: Product relations
-----------------------------

Particle-Hole channel (:math:`PH`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider the PH product

.. math::
   P(a\bar{b}\bar{c}d)
   & =
   \Gamma^{PH}(a\bar{b}u\bar{v}) \, \chi^{PH}_0(\bar{v}u\bar{c}d)
   \\
   & =
   \sum_{u\bar{v}}
   \iint_0^\beta d\tau_{u} d\tau_{\bar{v}} \,
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\tau_{a} \tau_{\bar{b}} \tau_{u} \tau_{\bar{v}})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d}(\tau_{\bar{v}} \tau_{u} \tau_{\bar{c}} \tau_{d}).

Fourier transforming :math:`\Gamma^{PH}` and :math:`\chi^{PH}_0` and explicitly inserting Kronecker delta functions for the total frequency conservation gives

.. math::
   P(a\bar{b}\bar{c}d)
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a + i \nu_{\bar{b}} \tau_{\bar{b}} + i \nu_{\bar{c}} \tau_{\bar{c}}- i \nu_{d} \tau_d
   \Big]
   \\ 
   & \times
   \frac{1}{\beta^2}
   \sum_{u \bar{v}}
   \sum_{\nu_{u} \nu_{\bar{v}}}
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\nu_a \nu_{\bar{b}} \nu_{u} \nu_{\bar{v}})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d}(\nu_{\bar{v}} \nu_u \nu_{\bar{c}} \nu_d)
   \\
   & \times
   \delta_{\nu_{a} - \nu_{\bar{b}} + \nu_{u} - \nu_{\bar{v}}, 0} 
   \delta_{\nu_{\bar{v}} - \nu_{u} + \nu_{\bar{c}} - \nu_{d}, 0}. 
   
Inserting the :math:`PH` frequency pairing of :eq:`matsubara_vertex` in this expression fulfills both Kronecker delta functions and reduce the summation by one frequency to

.. math::
   P(a\bar{b}\bar{c}d)
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu \tau_a + i (\nu + \omega) \tau_{\bar{b}} + i (\nu' + \omega) \tau_{\bar{c}} - i \nu' \tau_d
   \Big]
   \\
   & \times
   \frac{1}{\beta^2} \sum_{u \bar{v}} \sum_{\bar{\nu}}
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\nu, \nu+\omega, \bar{\nu} + \omega, \bar{\nu})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d}(\bar{\nu}, \bar{\nu} + \omega, \nu' + \omega, \nu').

Using the three frequency notation :math:`Q(\omega, \nu, \nu') \equiv Q(\nu, \nu+\omega, \nu'+\omega, \nu)` we get the final product relation

.. math::
   P^{PH}_{a\bar{b}\bar{c}d}(\omega, \nu,\nu')
   & =
   \frac{1}{\beta^2} \sum_{\bar{\nu} u\bar{v}}
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\omega,\nu, \bar{\nu})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d }(\omega,\bar{\nu}, \nu)
   \\
   & = 
   \frac{1}{\beta^2} \sum_{\bar{\nu} u\bar{v}}
   \Gamma^{PH}_{ \{ \nu, a\bar{b} \},\{ \bar{\nu}, \bar{v}u \}}(\omega)
   \,
   \chi^{PH}_{0, \{\bar{\nu}, \bar{v}u \},\{ \nu, d\bar{c} \}}(\omega).

.. note::

   The right hand side indices has to be permuted in order to make the product a direct matrix multiplication. I.e. the pairing reads

   .. math::
      P^{PH}_{abcd}(\omega, \nu, \nu') = P^{PH}_{\{\nu, ab \}, \{\nu', dc\}}(\omega).
   
..
   The orbital indices in the right term :math:`\chi^{PH}_0` are transposed, i.e, :math:`\{ \bar{\nu},\bar{v}u \}` and not :math:`\{ \bar{\nu}, u\bar{v} \}` as in our reference notes!
   This transpose has to be done in the index reordering when mapping to matrices!
   **I.e. the index ordering has to be DIFFERENT for the left and right hand side of the** :math:`PH` **product!**
   
Writing the reversed product :math:`P = \chi^{PH}_0 * \Gamma^{PH}` in slightly compressed notation we get

.. math::
   \mathcal{F} \big\{ P(\bar{a}bc\bar{d}) \big\} =
   \frac{1}{\beta^2} \sum_{\bar{u}v} \sum_{\bar{\nu}}
   \chi^{PH}_{0, \bar{a}b\bar{u}v}(\nu \nu+\omega, \bar{\nu} + \omega, \bar{\nu})
   \,
   \Gamma^{PH}_{v\bar{u}c\bar{d}}(\bar{\nu}, \bar{\nu} + \omega, \nu' + \omega, \nu').
   
where :math:`\mathcal{F}\{ \cdot \}` denotes Fourier transformation to four fermionic Matsubara frequency space. Thus, the product with grouped indices becomes

.. math::
   P_{\bar{a}bc\bar{d}}(\omega, \nu, \nu')
   & =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u}v}
   \chi^{PH}_{0, \bar{a}b\bar{u}v}(\omega, \nu, \bar{\nu})
   \,
   \Gamma^{PH}_{v\bar{u}c\bar{d}}(\omega, \bar{\nu}, \nu')
   \\
   & =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u}v}
   \chi^{PH}_{0, \{ \nu, \bar{a}b \}, \{\bar{\nu}, v\bar{u} \} }(\omega)
   \,
   \Gamma^{PH}_{\{ \bar{\nu} , v\bar{u} \}, \{ \nu', \bar{d}c\}}(\omega)
   
which shows that the same index grouping relations hold for both products :math:`\chi_0^{PH} * \Gamma^{PH}` and :math:`\Gamma^{PH} * \chi_0^{PH}`.

   
Vertical-Particle-Hole channel (:math:`\bar{PH}`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vertical-particle-hole product is defined in the channel's Bethe-Salpeter equation as

.. math::
   P(ab\bar{c}\bar{d})
   & =
   \Gamma^{\bar{PH}}(a\bar{u}v\bar{d})
   \,
   \chi_0^{\bar{PH}}(\bar{u}b\bar{c}v)
   \\
   & =
   \sum_{\bar{u}v} \iint_0^\beta d\tau_{\bar{u}} d\tau_v \,
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\tau_a, \tau_{\bar{u}}, \tau_v, \tau_{\bar{d}})
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\tau_{\bar{u}},\tau_b,\tau_{\bar{c}},\tau_v).

Fourier expansion gives

.. math::
   P(ab\bar{c}\bar{d})
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a + i \nu_{\bar{b}} \tau_{\bar{b}} + i \nu_{\bar{c}} \tau_{\bar{c}} - i \nu_{d} \tau_d
   \Big]
   \\
   & \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} v}
   \sum_{\nu_{\bar{u}} \nu_{v}}
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\nu_a \nu_{\bar{u}} \nu_v \nu_{\bar{d}})
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\nu_{\bar{u}} \nu_b \nu_{\bar{c}} \nu_v)
   \\
   & \times
   \delta_{\nu_a - \nu_{\bar{u}} + \nu_v - \nu_{\bar{d}}, 0}
   \delta_{\nu_{\bar{u}} - \nu_b + \nu_{\bar{c}} - \nu_v, 0}.

Inserting the :math:`\bar{PH}` channel frequency parametrization of Eq. :eq:`matsubara_vertex`, gives

.. math::
   P(ab\bar{c}\bar{d})
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu \tau_a + i \nu' \tau_{\bar{b}} + i (\nu' + \omega) \tau_{\bar{c}} - i (\nu + \omega) \tau_d
   \Big]
   \\
   & \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} v}
   \sum_{\bar{\nu}}
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\nu, \bar{\nu}, \bar{\nu} + \omega, \nu + \omega)
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\bar{\nu}, \nu', \nu' + \omega, \bar{\nu} + \omega)

using :math:`\bar{PH}` frequency notation and grouping indices we get

.. math::
   P_{ab\bar{c}\bar{d}}(\omega, \nu, \nu')
   & =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u} v}
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\omega, \nu, \bar{\nu})
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\omega, \bar{\nu}, \nu')
   \\
   & =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u} v}
   \Gamma^{\bar{PH}}_{\{ \nu, a\bar{d} \}, \{ \bar{\nu}, \bar{u}v \}}(\omega)
   \,
   \chi^{\bar{PH}}_{0, \{\bar{\nu}, \bar{u}v \}, \{\nu', b\bar{c} \} }(\omega).

The reversed product :math:`\chi^{\bar{PH}}_0 * \Gamma^{\bar{PH}}` can be analysed in the same way and gives the same index pairing.
   

Particle-Particle channel (:math:`PP`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
   P(abcd)
   & =
   \Gamma^{PP}(a\bar{u}c\bar{v}) 
   \,
   \chi^{PP}_0(\bar{u}b\bar{v}d)
   \\
   & =
   \sum_{\bar{u}\bar{v}}
   \iint_0^\beta d\tau_{\bar{u}} d\tau_{\bar{v}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}(\tau_a, \tau_{\bar{u}}, \tau_c, \tau_{\bar{v}}) 
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}(\tau_{\bar{u}}, \tau_b, \tau_{\bar{v}}, \tau_d).

Fourier transform
   
.. math::
   P(abcd)
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a - i \nu_{\bar{b}} \tau_{\bar{b}} - i \nu_{\bar{c}} \tau_{\bar{c}} - i \nu_{d} \tau_d
   \Big]
   \\
   & \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\nu_{\bar{u}} \nu_{\bar{v}}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}(\nu_a \nu_{\bar{u}} \nu_c \nu_{\bar{v}})
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}(\nu_{\bar{u}} \nu_b \nu_{\bar{v}} \nu_d)
   \\
   & \times
   \delta_{\nu_a - \nu_{\bar{u}} + \nu_c - \nu_{\bar{v}}, 0}
   \delta_{\nu_{\bar{u}} - \nu_b + \nu_{\bar{v}} - \nu_d, 0}.

Inserting the :math:`PP` channel frequency parametrization of Eq. :eq:`matsubara_vertex` gives
   
.. math::
   P(abcd)
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i(\nu) \tau_a - i (\nu') \tau_{\bar{b}}
   - i (\omega - \nu') \tau_{\bar{c}} - i (\omega - \nu') \tau_d
   \Big]
   \\
   & \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}
   (\nu, \bar{\nu}, \omega - \nu, \omega - \bar{\nu})
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}
   (\bar{\nu}, \nu', \omega - \bar{\nu}, \omega - \nu').

Collecting indices
   
.. math::
   P_{abcd}(\omega, \nu, \nu')
   & =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}
   (\omega, \nu, \bar{\nu})
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}
   (\omega, \bar{\nu}, \nu')
   \\
   & =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{ \{ \nu , ac \}, \{\bar{\nu}, \bar{u}\bar{v} \} }
   (\omega)
   \,
   \chi^{PP}_{0, \{ \bar{\nu}, \bar{u}\bar{v} \}, \{ \nu', bd\}}
   (\omega).
   
Crossed-Particle-Particle channel (:math:`PPx`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
   P(abcd)
   & =
   \Gamma^{PPx}(a\bar{u}c\bar{v}) 
   \
   \chi^{PPx}_0(\bar{v}b\bar{u}d)
   \\
   & =
   \sum_{\bar{u}\bar{v}}
   \iint_0^\beta d\tau_{\bar{u}} d\tau_{\bar{v}}
   \Gamma^{PPx}_{a\bar{u}c\bar{v}}(\tau_a, \tau_{\bar{u}}, \tau_c, \tau_{\bar{v}}) 
   \,
   \chi^{PPx}_{0, \bar{v}b\bar{u}d}(\tau_{\bar{v}}, \tau_b, \tau_{\bar{u}}, \tau_d).

Fourier transform
   
.. math::
   P(abcd)
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a - i \nu_{\bar{b}} \tau_{\bar{b}} - i \nu_{\bar{c}} \tau_{\bar{c}} - i \nu_{d} \tau_d
   \Big]
   \\
   & \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\nu_{\bar{u}} \nu_{\bar{v}}}
   \Gamma^{PPx}_{a\bar{u}c\bar{v}}(\nu_a \nu_{\bar{u}} \nu_c \nu_{\bar{v}})
   \,
   \chi^{PPx}_{0, \bar{v}b\bar{u}d}(\nu_{\bar{v}} \nu_b \nu_{\bar{u}} \nu_d)
   \\
   & \times
   \delta_{\nu_a - \nu_{\bar{u}} + \nu_c - \nu_{\bar{v}}, 0}
   \delta_{\nu_{\bar{v}} - \nu_b + \nu_{\bar{u}} - \nu_d, 0}.

Inserting the :math:`PPx` channel parametrization of Eq. :eq:`matsubara_vertex` gives

.. math::
   \nu_a - \nu_{\bar{u}} + \nu_c - \nu_{\bar{v}}
   & =
   \nu - \omega + \bar{\nu} + \omega - \nu - \bar{\nu}
   & = 0 \\
   \nu_{\bar{v}} - \nu_b + \nu_{\bar{u}} - \nu_d
   & =
   \bar{\nu} - \omega + \nu' + \omega - \bar{\nu} - \nu'
   & = 0,
   
.. math::
   P(abcd)
   & =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i(\nu) \tau_a - i (\omega - \nu') \tau_{\bar{b}}
   - i (\omega - \nu) \tau_{\bar{c}} - i (\nu') \tau_d
   \Big]
   \\
   & \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PPx}_{a\bar{u}c\bar{v}}
   (\nu, \omega - \bar{\nu}, \omega - \nu, \bar{\nu})
   \,
   \chi^{PPx}_{0, \bar{u}b\bar{v}d}
   (\bar{\nu}, \omega - \nu', \omega - \bar{\nu}, \nu').

Collecting indices
   
.. math::
   P^{PPx}_{abcd}(\omega, \nu, \nu')
   & =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}
   (\omega, \nu, \bar{\nu})
   \,
   \chi^{PP}_{0, \bar{v}b\bar{u}d}
   (\omega, \bar{\nu}, \nu')
   \\
   & =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PPx}_{ \{ \nu , ca \}, \{\bar{\nu}, \bar{u}\bar{v} \} }
   (\omega)
   \,
   \chi^{PPx}_{0, \{ \bar{\nu}, \bar{u}\bar{v} \}, \{ \nu', bd\}}
   (\omega).
   
.. note::

   The first index is permuted in the grouping, i.e.

   .. math::
      P^{PPx}_{abcd}(\omega, \nu, \nu')
      = P^{PPx}_{\{\nu, ca\}, \{ \nu', bd \}}(\omega).

