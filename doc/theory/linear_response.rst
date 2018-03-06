.. _linear_response:

Linear response from generalized susceptibility
===============================================

The static linear response can be obtained from the frequency dependent generalized susceptibility :math:`\chi^{(PH)}_{\bar{a}b\bar{c}d}` or the one-time dependent bosonic propagator :math:`\langle (\bar{a}b)(\tau) (\bar{c}d) \rangle` by summing/integrating over freqiency/time

.. math::
   \chi^{(PH)}_{\bar{a}b\bar{c}d}
   =
   \frac{1}{\beta} \int_0^\beta d\tau \, \langle (\bar{a}b)(\tau) (\bar{c}d) \rangle
   - \langle \bar{a} b \rangle \langle \bar{c} d \rangle
   =
   \frac{1}{\beta^2} \sum_{nm}
   \chi^{(PH)}_{\bar{a}b\bar{c}d}(\omega_0, \nu_n, \nu_m)
   

Linear response from applied external field
===========================================

The linear response of a system is related to the systems response to the application of an external field. Hence, a brute force approach to compute the static linear response is to apply a general quadratic field :math:`F_{\bar{a}b} ( \bar{a} b + \bar{b} a)` to the system and compute the derivative of the single particle density matrix :math:`\rho_{\bar{a}b} \equiv \langle \bar{a}b \rangle` with respect to the field :math:`F_{\bar{a}b}` in the limit :math:`F_{\bar{a}b} \rightarrow 0`. In other words we can compute the system response :math:`R_{\bar{a}b\bar{c}d}` as

.. math::
   R_{\bar{a}b\bar{c}d} \equiv
   \partial_{F_{\bar{a}b}} \langle \bar{c} d \rangle |_{F_{\bar{a}b} = 0}

However, the central difference between the system response :math:`R_{\bar{a}b\bar{c}d}` and the static the generalized susceptibility :math:`\chi^{(PH)}_{\bar{a}b\bar{c}d}` is that the applied field has to be Hermitian. Whence, they are related by

.. math::
   R_{\bar{a}b\bar{c}d} = \chi^{(PH)}_{\bar{a}b\bar{c}d} + \chi^{(PH)}_{\bar{b}a\bar{c}d}

and the response :math:`R_{\bar{a}b\bar{c}d}` is symmetric in :math:`\bar{a}b` transposes :math:`R_{\bar{a}b\bar{c}d} = R_{\bar{b}a\bar{c}d}`.
   

Reconstructing :math:`\chi` from :math:`R`
------------------------------------------

from the symmetry relation

.. math::
   \langle (\bar{a}b)(\tau) (\bar{c}d) \rangle
   =
   [ \xi \langle (\bar{d}c)(\beta - \tau) (\bar{b}a) \rangle ]^*

we can compute four relations involving pair-wise permutations of :math:`\chi^{(PH)}_{\bar{a}b\bar{c}d}` from :math:`R_{\bar{a}b\bar{c}d}`

.. math::
   R_{\bar{a}b\bar{c}d} = \chi^{(PH)}_{\bar{a}b\bar{c}d} + \chi^{(PH)}_{\bar{b}a\bar{c}d} \\
   R_{\bar{a}b\bar{d}c} = \chi^{(PH)}_{\bar{a}b\bar{d}c} + \chi^{(PH)}_{\bar{b}a\bar{d}c} \\
   [R_{\bar{c}d\bar{a}b}]^* = \chi^{(PH)}_{\bar{b}a\bar{d}c} + \chi^{(PH)}_{\bar{b}a\bar{c}d} \\
   [R_{\bar{c}d\bar{b}a}]^* = \chi^{(PH)}_{\bar{a}b\bar{d}c} + \chi^{(PH)}_{\bar{a}b\bar{c}d} \\

However, these relations are not linearly independent and can not be solved for the :math:`\chi^{(PH)}` components.
