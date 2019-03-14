.. _derivation_index_pairing:

Derivation: Product relations
=============================

.. _derivation_index_pairing_ph:

Particle-Hole channel (:math:`PH`)
----------------------------------

Consider the PH product

.. math::
   \begin{multline}
   P(a\bar{b}\bar{c}d) =
   \Gamma^{PH}(a\bar{b}u\bar{v}) \, \chi^{PH}_0(\bar{v}u\bar{c}d)
   \\ =
   \sum_{u\bar{v}}
   \iint_0^\beta d\tau_{u} d\tau_{\bar{v}} \,
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\tau_{a} \tau_{\bar{b}} \tau_{u} \tau_{\bar{v}})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d}(\tau_{\bar{v}} \tau_{u} \tau_{\bar{c}} \tau_{d})
   \end{multline}

Fourier transforming :math:`\Gamma^{PH}` and :math:`\chi^{PH}_0` and explicitly inserting Kronecker delta functions for the total frequency conservation gives

.. math::
   P(a\bar{b}\bar{c}d) =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a + i \nu_{\bar{b}} \tau_{\bar{b}} + i \nu_{\bar{c}} \tau_{\bar{c}}- i \nu_{d} \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2}
   \sum_{u \bar{v}}
   \sum_{\nu_{u} \nu_{\bar{v}}}
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\nu_a \nu_{\bar{b}} \nu_{u} \nu_{\bar{v}})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d}(\nu_{\bar{v}} \nu_u \nu_{\bar{c}} \nu_d)
   \\ \times
   \delta_{\nu_{a} - \nu_{\bar{b}} + \nu_{u} - \nu_{\bar{v}}, 0} 
   \delta_{\nu_{\bar{v}} - \nu_{u} + \nu_{\bar{c}} - \nu_{d}, 0} 
   
Inserting the :math:`PH` frequency pairing :eq:`PH_freq` in this expression fulfills both Kronecker delta functions and reduce the summation by one frequency to

.. math::
   P(a\bar{b}\bar{c}d) =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu \tau_a + i (\nu + \omega) \tau_{\bar{b}} + i (\nu' + \omega) \tau_{\bar{c}} - i \nu' \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2} \sum_{u \bar{v}} \sum_{\bar{\nu}}
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\nu, \nu+\omega, \bar{\nu} + \omega, \bar{\nu})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d}(\bar{\nu}, \bar{\nu} + \omega, \nu' + \omega, \nu')

Using the three frequency notation :math:`Q(\omega, \nu, \nu') \equiv Q(\nu, \nu+\omega, \nu'+\omega, \nu)` we get the final product relation

.. math::
   P^{PH}_{a\bar{b}\bar{c}d}(\omega, \nu,\nu') =
   \frac{1}{\beta^2} \sum_{\bar{\nu} u\bar{v}}
   \Gamma^{PH}_{a\bar{b}u\bar{v}}(\omega,\nu, \bar{\nu})
   \,
   \chi^{PH}_{0, \bar{v}u\bar{c}d }(\omega,\bar{\nu}, \nu)
   \\ = 
   \frac{1}{\beta^2} \sum_{\bar{\nu} u\bar{v}}
   \Gamma^{PH}_{ \{ \nu, a\bar{b} \},\{ \bar{\nu}, \bar{v}u \}}(\omega)
   \,
   \chi^{PH}_{0, \{\bar{\nu}, \bar{v}u \},\{ \nu, d\bar{c} \}}(\omega)

.. note::

   The right hand side indices has to be permuted in order to make the product a direct matrix multiplication. I.e. the pairing reads

   .. math::
      P^{PH}_{abcd}(\omega, \nu, \nu') = P^{PH}_{\{\nu, ab \}, \{\nu', dc\}}(\omega)
   
..
   The orbital indices in the right term :math:`\chi^{PH}_0` are transposed, i.e, :math:`\{ \bar{\nu},\bar{v}u \}` and not :math:`\{ \bar{\nu}, u\bar{v} \}` as in our reference notes!

   This transpose has to be done in the index reordering when mapping to matrices!

   **I.e. the index ordering has to be DIFFERENT for the left and right hand side of the** :math:`PH` **product!**
   
Writing the reversed product :math:`P = \chi^{PH}_0 * \Gamma^{PH}` in slightly compressed notation we get

.. math::
   \mathcal{F} \big\{ P(\bar{a}bc\bar{d}) \big\}
   \\ =
   \frac{1}{\beta^2} \sum_{\bar{u}v} \sum_{\bar{\nu}}
   \chi^{PH}_{0, \bar{a}b\bar{u}v}(\nu \nu+\omega, \bar{\nu} + \omega, \bar{\nu})
   \,
   \Gamma^{PH}_{v\bar{u}c\bar{d}}(\bar{\nu}, \bar{\nu} + \omega, \nu' + \omega, \nu')
   
where :math:`\mathcal{F}\{ \cdot \}` denotes Fourier transformation to four fermionic Matsubara frequency space. Thus, the product with grouped indices becomes

.. math::
   P_{\bar{a}bc\bar{d}}(\omega, \nu, \nu')
   =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u}v}
   \chi^{PH}_{0, \bar{a}b\bar{u}v}(\omega, \nu, \bar{\nu})
   \,
   \Gamma^{PH}_{v\bar{u}c\bar{d}}(\omega, \bar{\nu}, \nu')
   \\=
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u}v}
   \chi^{PH}_{0, \{ \nu, \bar{a}b \}, \{\bar{\nu}, v\bar{u} \} }(\omega)
   \,
   \Gamma^{PH}_{\{ \bar{\nu} , v\bar{u} \}, \{ \nu', \bar{d}c\}}(\omega)
   
which shows that the same index grouping relations hold for both products :math:`\chi_0^{PH} * \Gamma^{PH}` and :math:`\Gamma^{PH} * \chi_0^{PH}`.

   
Vertical-Particle-Hole channel (:math:`\bar{PH}`)
-------------------------------------------------

The vertical-particle-hole product is defined in the channel's Bethe-Salpeter equation as

.. math::
   \begin{multline}
   P(ab\bar{c}\bar{d}) =
   \Gamma^{\bar{PH}}(a\bar{u}v\bar{d})
   \,
   \chi_0^{\bar{PH}}(\bar{u}b\bar{c}v)
   \\ =
   \sum_{\bar{u}v} \iint_0^\beta d\tau_{\bar{u}} d\tau_v \,
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\tau_a, \tau_{\bar{u}}, \tau_v, \tau_{\bar{d}})
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\tau_{\bar{u}},\tau_b,\tau_{\bar{c}},\tau_v)
   \end{multline}

Fourier expansion gives

.. math::
   P(ab\bar{c}\bar{d}) =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a + i \nu_{\bar{b}} \tau_{\bar{b}} + i \nu_{\bar{c}} \tau_{\bar{c}} - i \nu_{d} \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} v}
   \sum_{\nu_{\bar{u}} \nu_{v}}
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\nu_a \nu_{\bar{u}} \nu_v \nu_{\bar{d}})
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\nu_{\bar{u}} \nu_b \nu_{\bar{c}} \nu_v)
   \\ \times
   \delta_{\nu_a - \nu_{\bar{u}} + \nu_v - \nu_{\bar{d}}, 0}
   \delta_{\nu_{\bar{u}} - \nu_b + \nu_{\bar{c}} - \nu_v, 0}

Inserting the :math:`\bar{PH}` channel frequency parametrization of Eq. :eq:`PHbar_freq`, gives

.. math::
   P(ab\bar{c}\bar{d}) =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu \tau_a + i \nu' \tau_{\bar{b}} + i (\nu' + \omega) \tau_{\bar{c}} - i (\nu + \omega) \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} v}
   \sum_{\bar{\nu}}
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\nu, \bar{\nu}, \bar{\nu} + \omega, \nu + \omega)
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\bar{\nu}, \nu', \nu' + \omega, \bar{\nu} + \omega)

using :math:`\bar{PH}` frequency notation and grouping indices we get

.. math::
   P_{ab\bar{c}\bar{d}}(\omega, \nu, \nu')
   =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u} v}
   \Gamma^{\bar{PH}}_{a\bar{u}v\bar{d}}(\omega, \nu, \bar{\nu})
   \,
   \chi^{\bar{PH}}_{0, \bar{u}b\bar{c}v}(\omega, \bar{\nu}, \nu')
   \\ =
   \frac{1}{\beta^2} \sum_{\bar{\nu}, \bar{u} v}
   \Gamma^{\bar{PH}}_{\{ \nu, a\bar{d} \}, \{ \bar{\nu}, \bar{u}v \}}(\omega)
   \,
   \chi^{\bar{PH}}_{0, \{\bar{\nu}, \bar{u}v \}, \{\nu', b\bar{c} \} }(\omega)

The reversed product :math:`\chi^{\bar{PH}}_0 * \Gamma^{\bar{PH}}` can be analysed in the same way and gives the same index pairing.
   

Particle-Particle channel (:math:`PP`)
--------------------------------------

.. math::
   \begin{multline}
   P(abcd) =
   \Gamma^{PP}(a\bar{u}c\bar{v}) 
   \,
   \chi^{PP}_0(\bar{u}b\bar{v}d)
   \\ =
   \sum_{\bar{u}\bar{v}}
   \iint_0^\beta d\tau_{\bar{u}} d\tau_{\bar{v}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}(\tau_a, \tau_{\bar{u}}, \tau_c, \tau_{\bar{v}}) 
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}(\tau_{\bar{u}}, \tau_b, \tau_{\bar{v}}, \tau_d)
   \end{multline}

Fourier transform
   
.. math::
   P(abcd)
   =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a - i \nu_{\bar{b}} \tau_{\bar{b}} - i \nu_{\bar{c}} \tau_{\bar{c}} - i \nu_{d} \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\nu_{\bar{u}} \nu_{\bar{v}}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}(\nu_a \nu_{\bar{u}} \nu_c \nu_{\bar{v}})
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}(\nu_{\bar{u}} \nu_b \nu_{\bar{v}} \nu_d)
   \\ \times
   \delta_{\nu_a - \nu_{\bar{u}} + \nu_c - \nu_{\bar{v}}, 0}
   \delta_{\nu_{\bar{u}} - \nu_b + \nu_{\bar{v}} - \nu_d, 0}

Inserting Eq. :eq:`PP_freq` gives
   
.. math::
   P(abcd)
   =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i(\nu) \tau_a - i (\nu') \tau_{\bar{b}}
   - i (\omega - \nu') \tau_{\bar{c}} - i (\omega - \nu') \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}
   (\nu, \bar{\nu}, \omega - \nu, \omega - \bar{\nu})
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}
   (\bar{\nu}, \nu', \omega - \bar{\nu}, \omega - \nu')

Collecting indices
   
.. math::
   P_{abcd}(\omega, \nu, \nu')
   =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}
   (\omega, \nu, \bar{\nu})
   \,
   \chi^{PP}_{0, \bar{u}b\bar{v}d}
   (\omega, \bar{\nu}, \nu')
   \\ =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{ \{ \nu , ac \}, \{\bar{\nu}, \bar{u}\bar{v} \} }
   (\omega)
   \,
   \chi^{PP}_{0, \{ \bar{\nu}, \bar{u}\bar{v} \}, \{ \nu', bd\}}
   (\omega)
   
Crossed-Particle-Particle channel (:math:`PPx`)
-----------------------------------------------

.. math::
   \begin{multline}
   P(abcd) =
   \Gamma^{PPx}(a\bar{u}c\bar{v}) 
   \
   \chi^{PPx}_0(\bar{v}b\bar{u}d)
   \\ =
   \sum_{\bar{u}\bar{v}}
   \iint_0^\beta d\tau_{\bar{u}} d\tau_{\bar{v}}
   \Gamma^{PPx}_{a\bar{u}c\bar{v}}(\tau_a, \tau_{\bar{u}}, \tau_c, \tau_{\bar{v}}) 
   \,
   \chi^{PPx}_{0, \bar{v}b\bar{u}d}(\tau_{\bar{v}}, \tau_b, \tau_{\bar{u}}, \tau_d)
   \end{multline}

Fourier transform
   
.. math::
   P(abcd)
   =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i\nu_a \tau_a - i \nu_{\bar{b}} \tau_{\bar{b}} - i \nu_{\bar{c}} \tau_{\bar{c}} - i \nu_{d} \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\nu_{\bar{u}} \nu_{\bar{v}}}
   \Gamma^{PPx}_{a\bar{u}c\bar{v}}(\nu_a \nu_{\bar{u}} \nu_c \nu_{\bar{v}})
   \,
   \chi^{PPx}_{0, \bar{v}b\bar{u}d}(\nu_{\bar{v}} \nu_b \nu_{\bar{u}} \nu_d)
   \\ \times
   \delta_{\nu_a - \nu_{\bar{u}} + \nu_c - \nu_{\bar{v}}, 0}
   \delta_{\nu_{\bar{v}} - \nu_b + \nu_{\bar{u}} - \nu_d, 0}

Inserting Eq. :eq:`PPx_freq` gives

.. math::
   \nu_a - \nu_{\bar{u}} + \nu_c - \nu_{\bar{v}}
   =
   \nu - \omega + \bar{\nu} + \omega - \nu - \bar{\nu}
   = 0 \\
   \nu_{\bar{v}} - \nu_b + \nu_{\bar{u}} - \nu_d
   =
   \bar{\nu} - \omega + \nu' + \omega - \bar{\nu} - \nu' = 0
   
.. math::
   P(abcd)
   =
   \frac{1}{\beta^4} \sum
   \exp \Big[
   -i(\nu) \tau_a - i (\omega - \nu') \tau_{\bar{b}}
   - i (\omega - \nu) \tau_{\bar{c}} - i (\nu') \tau_d
   \Big]
   \\ \times
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PPx}_{a\bar{u}c\bar{v}}
   (\nu, \omega - \bar{\nu}, \omega - \nu, \bar{\nu})
   \,
   \chi^{PPx}_{0, \bar{u}b\bar{v}d}
   (\bar{\nu}, \omega - \nu', \omega - \bar{\nu}, \nu')

Collecting indices
   
.. math::
   P^{PPx}_{abcd}(\omega, \nu, \nu')
   =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PP}_{a\bar{u}c\bar{v}}
   (\omega, \nu, \bar{\nu})
   \,
   \chi^{PP}_{0, \bar{v}b\bar{u}d}
   (\omega, \bar{\nu}, \nu')
   \\ =
   \frac{1}{\beta^2}
   \sum_{\bar{u} \bar{v}}
   \sum_{\bar{\nu}}
   \Gamma^{PPx}_{ \{ \nu , ca \}, \{\bar{\nu}, \bar{u}\bar{v} \} }
   (\omega)
   \,
   \chi^{PPx}_{0, \{ \bar{\nu}, \bar{u}\bar{v} \}, \{ \nu', bd\}}
   (\omega)
   
.. note::

   The first index is permuted in the grouping, i.e.

   .. math::
      P^{PPx}_{abcd}(\omega, \nu, \nu')
      = P^{PPx}_{\{\nu, ca\}, \{ \nu', bd \}}(\omega)
