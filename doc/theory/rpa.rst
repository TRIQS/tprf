.. _rpa:

Random phase approximation (RPA)
================================

Interaction
-----------

.. math::
   H_{int} =
   \sum_{abcd} V(\bar{a}\bar{b}cd) \,\, \bar{a} \bar{b} c d

Hermicity
   
.. math::
   H_{int} = H_{int}^{\dagger}

.. math::
   \Big[ V(\bar{a}\bar{b}cd) \,\, \bar{a}\bar{b}cd  \Big]^\dagger
   =
   V^*(\bar{a}\bar{b}cd) \,\, \bar{d} \bar{c} b a
   =
   V(\bar{d}\bar{c}ba) \,\, \bar{d} \bar{c} b a
   
.. math::
   V^*(\bar{a}\bar{b}cd) = V(\bar{d}\bar{c}ba)

.. math::
   V(\bar{a}\bar{b}cd) = V^*(\bar{d}\bar{c}ba)

.. note::

   **Techincal detail:** In TRIQS :math:`H_{int}` is represented with normal ordered operators where the indices of the operators are sorted in lexicographical order for the creation operators and reverse lexicographical order for the annihilation operators.

   For this reason the interaction tensor :math:`\tilde{V}(\bar{a}\bar{b}cd)` obtained from :math:`H_{int}` has to be "up-folded" to all permutations of annihilation and creation operator pairs, respectively

   .. math::
      V(\bar{a}\bar{b}cd) =
      \frac{1}{4}
      \Big[
        \tilde{V}(\bar{a}\bar{b}cd)
        - \tilde{V}(\bar{b}\bar{a}cd)
	+ \tilde{V}(\bar{b}\bar{a}dc)
	- \tilde{V}(\bar{a}\bar{b}dc)
      \Big]

RPA susceptibility in the particle-hole channel (PH)
----------------------------------------------------

.. math::
   \chi^{(PH)}(\bar{a}b\bar{c}d) =
   \chi_0(\bar{a}b\bar{c}d)
   + \chi_0(\bar{a}b\bar{p}q) \,
     \Gamma^{(PH)}(q\bar{p}s\bar{r}) \,
     \chi^{(PH)}(\bar{r}s\bar{c}d)

RPA approximation in the particle hole channel (PH)

.. math::
   \Gamma^{(PH)}(q\bar{p}s\bar{r}) \approx \pm V(\bar{p}\bar{r}qs)

.. note::
   Figure out the sign in the prefactor of the RPA approximation for :math:`\Gamma^{(PH)}`.
   
.. math::
   \chi^{(PH)}_{RPA}(\bar{a}b\bar{c}d) =
   \chi_0(\bar{a}b\bar{c}d)
   + \chi_0(\bar{a}b\bar{p}q) \,
     V(\bar{p}\bar{r}qs) \,
     \chi^{(PH)}_{RPA}(\bar{r}s\bar{c}d)

**Alternative 1**

.. math::
   \chi = \chi_0 + \chi_0 V \chi
     
.. math::
   \chi_0^{-1} = \chi^{-1} + V

.. math::
   \chi_0^{-1} - V = \chi^{-1} 

.. math::
   [1 - V \chi_0] \, \chi_0^{-1} = \chi^{-1} 

.. math::
   \chi_0 \, [1 - V \chi_0]^{-1} = \chi


**Alternative start** (becomes equivalent in 2nd step, ok)

.. math::
   \chi = \chi_0 + \chi V \chi_0

.. math::
   \chi_0^{-1} = \chi^{-1} + V


**Alternative ending** 

.. math::
   \chi_0^{-1} - V = \chi^{-1} 

.. math::
   \chi_0^{-1} \, [ 1 - \chi_0 V ] = \chi^{-1}
   
.. math::
   [ 1 - \chi_0 V ]^{-1} \, \chi_0 = \chi
   

\Tr[ V ... V ]
