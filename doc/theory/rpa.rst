.. _rpa:

Random phase approximation (RPA)
================================

Interaction
-----------

In `triqs` the interaction Hamiltonian is represented as a sum of monomials of quartic operators

.. math::
   H_{int} =
   \sum_{ \{\bar{a}\bar{b}cd \}_s} V(\bar{a}\bar{b}cd) \,\, \bar{a} \bar{b} c d

where the sum runs over all unique sets of :math:`\bar{a}\bar{b}cd` of normal ordered and lexicographically ordered operators.

.. note::

   This is a unique representation of the Hamiltonian and a unique representation of the prefactor :math:`V(\bar{a}\bar{b}cd)`, in contrast to representations where we allow any permutation of :math:`\bar{a}\bar{b}` and :math:`cd`.

RPA tensor
----------

In RPA we approximate the vertex :math:`\Gamma` in the Bethe-Salpeter equation

.. math::
   \chi^{(PH)}(\bar{a}b\bar{c}d) =
   \chi_0(\bar{a}b\bar{c}d)
   + \chi_0(\bar{a}b\bar{p}q) \,
     \Gamma^{(PH)}(q\bar{p}s\bar{r}) \,
     \chi^{(PH)}(\bar{r}s\bar{c}d)

by a constant rank 4 tensor :math:`U(\bar{a}b\bar{c}d)`
     
.. math::
   \Gamma^{(PH)}(q\bar{p}s\bar{r}) \approx U(q\bar{p}s\bar{r})

To determine the relation between :math:`U(\bar{a}b\bar{c}d)` and :math:`V(\bar{a}\bar{b}cd)` we expand :math:`\chi` to first order in :math:`V`

The generalized susceptibility is defined as

.. math::
   \chi(\bar{a}b\bar{c}d) =
   \langle \bar{a}b\bar{c}d \rangle
   - \langle b \bar{a} \rangle \langle d \bar{c} \rangle

to zeroth order in :math:`V` we get the bare susceptibility
     
.. math::
   [\chi]_0 = \chi_0 = - \langle d\bar{a} \rangle \langle b \bar{c} \rangle

the first order is given by
   
.. math::
   [\chi]_1 =
   - \langle \bar{a}b\bar{c}d H_{int} \rangle
   + \langle b \bar{a} H_{int} \rangle \langle d \bar{c} \rangle
   + \langle b \bar{a} \rangle \langle d \bar{c} H_{int} \rangle
   \\ =
   \sum_{ \{\bar{A}\bar{B}CD \}_s }
   V(\bar{A}\bar{B}CD)
   \langle \bar{a}\bar{c} CD \rangle \langle bd \bar{A}\bar{B} \rangle

where we in the last step perform the restricted summation over unique interaction terms, as defined above, and use the fact that all contractions of :math:`d\bar{c}` and :math:`b\bar{a}` in the first term are canceled by the two last terms.

Performing the Wick contraction of the result and pairing the quadratic expectation values into :math:`\chi_0` terms gives

.. math::
   [\chi]_1 =
   \sum_{ \{ \bar{A}\bar{B}CD \}_s}
   V(\bar{A}\bar{B}CD)
   \Big[
     \chi_0(\bar{a}b \bar{A}C) \chi_0(\bar{B}D\bar{c}d) 
   - \chi_0(\bar{a}b \bar{A}D) \chi_0(\bar{B}C\bar{c}d) \\
   - \chi_0(\bar{a}b \bar{B}C) \chi_0(\bar{A}D\bar{c}d)
   + \chi_0(\bar{a}b \bar{B}D) \chi_0(\bar{A}C\bar{c}d)
   \Big]

by defining the tensor :math:`U` as

.. math::
   U(\bar{A}C\bar{B}D) = + V(\bar{A}\bar{B}CD)\\
   U(\bar{A}D\bar{B}C) = - V(\bar{A}\bar{B}CD)\\
   U(\bar{B}C\bar{A}D) = - V(\bar{A}\bar{B}CD)\\
   U(\bar{B}D\bar{A}C) = + V(\bar{A}\bar{B}CD)

we can rewrite the above equation as an unrestricted sum over :math:`U(\bar{A}B\bar{C}D)`
   
.. math::
   [\chi]_1 =
   \sum_{ \bar{A}B\bar{C}D }
     \chi_0(\bar{a}b\bar{A}B)
     U(\bar{A}B\bar{C}D)
     \chi_0(\bar{C}D\bar{c}d)

which determines that the RPA :math:`U(\bar{A}C\bar{B}D)` tensor transforms as the prefactor of

.. math::
   -V(\bar{A}\bar{B}CD) \bar{A}\bar{B}CD

under permutations of the indices.
