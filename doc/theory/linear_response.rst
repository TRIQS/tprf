.. _linear_response:

Linear response from applied external field
===========================================

The linear response of a system is related to the systems response to the application of an external field. Hence, a brute force approach to compute the static linear response is to apply a general quadratic field :math:`F_{\bar{a}b} ( \bar{a} b + \bar{b} a)` to the system and compute the derivative of the single particle density matrix with respect to the field :math:`F_{\bar{a}b}` in the limit :math:`F_{\bar{a}b} \rightarrow 0`. In other words we can compute the system response :math:`R_{\bar{a}b\bar{c}d}}` as

.. math::
   R_{\bar{a}b\bar{c}d}} \equiv
   \partial_{F_{\bar{a}b}} \langle \bar{c} d \rangle |_{F_{\bar{a}b} = 0}
