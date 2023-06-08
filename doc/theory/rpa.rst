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

In RPA we approximate the vertex :math:`\Gamma^{PH}` in the Bethe-Salpeter equation

.. math::
   \chi^{PH}(\bar{a}b\bar{c}d) =
   \chi_0^{PH}(\bar{a}b\bar{c}d)
   + \chi_0^{PH}(\bar{a}b\bar{p}q) \,
     \Gamma^{PH}(q\bar{p}s\bar{r}) \,
     \chi^{PH}(\bar{r}s\bar{c}d)

by a constant rank 4 tensor :math:`U(a\bar{b}c\bar{d})`
     
.. math::
   \Gamma^{PH}(q\bar{p}s\bar{r}) \approx U(q\bar{p}s\bar{r})

To determine the relation between :math:`U(a\bar{b}c\bar{d})` and :math:`V(\bar{a}\bar{b}cd)` we expand :math:`\chi` to first order in :math:`V`.
The generalized susceptibility is defined as

.. math::
   \chi(\bar{a}b\bar{c}d) =
   \langle \bar{a}b\bar{c}d \rangle
   - \langle b \bar{a} \rangle \langle d \bar{c} \rangle.

To zeroth order in :math:`V`, we get the bare susceptibility
     
.. math::
   [\chi]_0 = \chi_0 = - \langle d\bar{a} \rangle \langle b \bar{c} \rangle.

The first order is given by
   
.. math::
   [\chi]_1 =
   - \langle \bar{a}b\bar{c}d H_{int} \rangle
   + \langle b \bar{a} H_{int} \rangle \langle d \bar{c} \rangle
   + \langle b \bar{a} \rangle \langle d \bar{c} H_{int} \rangle
   \\ =
   \sum_{ \{\bar{A}\bar{B}CD \}_s }
   V(\bar{A}\bar{B}CD)
   \langle \bar{a}\bar{c} CD \rangle \langle bd \bar{A}\bar{B} \rangle

where in the last step we performed the restricted summation over unique interaction terms, as defined above, and use the fact that all contractions of :math:`d\bar{c}` and :math:`b\bar{a}` in the first term are canceled by the two last terms.

Performing the Wick contraction and pairing the quadratic expectation values into :math:`\chi_0` terms gives

.. math::
   [\chi]_1 =
   \sum_{ \{ \bar{A}\bar{B}CD \}_s}
   V(\bar{A}\bar{B}CD)
   \Big[
     \chi_0(\bar{a}b \bar{A}C) \chi_0(\bar{B}D\bar{c}d) 
   - \chi_0(\bar{a}b \bar{A}D) \chi_0(\bar{B}C\bar{c}d) \\
   - \chi_0(\bar{a}b \bar{B}C) \chi_0(\bar{A}D\bar{c}d)
   + \chi_0(\bar{a}b \bar{B}D) \chi_0(\bar{A}C\bar{c}d)
   \Big].

Defining the tensor :math:`U` as

.. math::
   U(\bar{A}C\bar{B}D) = + V(\bar{A}\bar{B}CD)\\
   U(\bar{A}D\bar{B}C) = - V(\bar{A}\bar{B}CD)\\
   U(\bar{B}C\bar{A}D) = - V(\bar{A}\bar{B}CD)\\
   U(\bar{B}D\bar{A}C) = + V(\bar{A}\bar{B}CD),

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


.. _matrix_rpa:

Spin-independent RPA calculations
---------------------------------

In many papers, such as `Takimoto, et. al., PRB 69, 104504 (2004)
<https://arxiv.org/abs/cond-mat/0309575>`_, people use a matrix scheme to calculate 
spin- and charge-susceptibilites in RPA.
In contrast to the generalized susceptibility used in TPRF they don't use an explicit
spin index, which is valid in SU(2) symmetric systems.
This allows for more efficient calculations both in memory and time.
A mapping between both methods is therefore advantageous and shown in
this documentation entry.

With the mapping we can obtain a full generalized RPA susceptibility from
spin-independent calculations by decoupling a spin-dependent interaction
into spin and charge channel.

Matrix RPA
^^^^^^^^^^

The matrix RPA is based on two equations,
one for the spin-susceptibility

.. math::
    \hat{\chi}^{(\mathrm{s})} = 
    \big(\mathbb{1} - \hat{\chi}^{(0)} \hat{U}^{(\mathrm{s})}\big)^{-1}
    \hat{\chi}^{(0)}\,,

and the other for the charge susceptibility

.. math::
    \hat{\chi}^{(\mathrm{c})} = 
    \big(\mathbb{1} + \hat{\chi}^{(0)} \hat{U}^{(\mathrm{c})}\big)^{-1}
    \hat{\chi}^{(0)}\,.

.. note::
    This is different to the equations in the paper `Takimoto, et. al., PRB 69, 104504 (2004)
    <https://arxiv.org/abs/cond-mat/0309575>`_ where they have :math:`\chi` and :math:`U` flipped.
    But I think that this is the correct way to write it.
    This difference must be resolved, because it will have difference on the particle-particle
    vertex, and therefore calculations of the Eliashberg equation in the RPA limit.

Here every quantity with a hat is a matrix and :math:`\hat{\chi}^{(0)}` is therefore
the matrix representation of the bare particle-hole bubble and the same holds for the
spin channel :math:`\hat{U}^{(\mathrm{s})}` and charge channel :math:`\hat{U}^{(\mathrm{c})}`
of the interaction.

The matrices are given as

.. math::
    \hat{\chi}^{(0)} = 
    \begin{pmatrix}
    \chi^{(0)}_{0000} & \chi^{(0)}_{0011} & \chi^{(0)}_{0001} & \chi^{(0)}_{0010}\\
    \chi^{(0)}_{1100} & \chi^{(0)}_{1111} & \chi^{(0)}_{1101} & \chi^{(0)}_{1110}\\
    \chi^{(0)}_{0100} & \chi^{(0)}_{0111} & \chi^{(0)}_{0101} & \chi^{(0)}_{0110}\\
    \chi^{(0)}_{1000} & \chi^{(0)}_{1011} & \chi^{(0)}_{1001} & \chi^{(0)}_{1010}\\
    \end{pmatrix}\,,

and for a Kanamori interaction by

.. math::
    \hat{U}^{(\mathrm{s})} = 
    \begin{pmatrix}
    U & J & 0 & 0\\
    J & U & 0 & 0\\
    0 & 0 & U' & J'\\
    0 & 0 & J' & U'\\
    \end{pmatrix}\,,

.. math::
    \hat{U}^{(\mathrm{c})} = 
    \begin{pmatrix}
    U & 2U'-J & 0 & 0\\
    2U'-J & U & 0 & 0\\
    0 & 0 & -U'+2J & J'\\
    0 & 0 & J' & -U'+2J\\
    \end{pmatrix}\,.

The equation for the elements of the bare-particle hole bubbles is given by

.. math::
   \chi^{(0)}_{\bar{\alpha}\beta\gamma\bar{\delta}}(q, i\omega) = -\frac{T}{N} \sum_{k, \nu} 
   G_{\gamma\bar{\alpha}}(k+q, i\nu + i\omega) G_{\beta\bar{\delta}}(k, i\nu)  \,,

were we used greek indices, which we will do exclusively for matrix RPA quantities to
highlight that they are spin-independent.
The spin and charge channel interaction are given by

.. math::
    U^{(\mathrm{s})}(\alpha\bar{\beta}\bar{\gamma}\delta) =
    \begin{cases}
    U, & \mathrm{if}\;\alpha=\bar{\beta}=\bar{\gamma}=\delta \\
    U', & \mathrm{if}\;\alpha=\bar{\gamma}\neq \bar{\beta}=\delta \\
    J, & \mathrm{if}\;\alpha=\bar{\beta}\neq \bar{\gamma}=\delta \\
    J', & \mathrm{if}\;\alpha=\delta\neq \bar{\beta}=\bar{\gamma} \\
    0, & \mathrm{else}
    \end{cases}\,,

.. math::
    U^{(\mathrm{c})}(\alpha\bar{\beta}\bar{\gamma}\delta) =
    \begin{cases}
    U, & \mathrm{if}\;\alpha=\bar{\beta}=\bar{\gamma}=\delta \\
    -U'+2J, & \mathrm{if}\;\alpha=\bar{\gamma}\neq \bar{\beta}=\delta \\
    2U'-J, & \mathrm{if}\;\alpha=\bar{\beta}\neq \bar{\gamma}=\delta \\
    J', & \mathrm{if}\;\alpha=\delta\neq \bar{\beta}=\bar{\gamma} \\
    0, & \mathrm{else}
    \end{cases}\,.

Here we have to note that operator order used in the matrix RPA differs from the one we use in TPRF.
While in TPRF the susceptibilites are ordered as :math:`c^\dagger cc^\dagger c` and the vertices
as :math:`cc^\dagger cc^\dagger`, in matrix RPA the last two indices are flipped.
This means :math:`c^\dagger ccc^\dagger` for susceptibilites and :math:`cc^\dagger c^\dagger c`
for vertices.
This flipping of the last two indices corresponds to the particle-hole product, see
:ref:`derivation_index_pairing`, which the matrix RPA explicitly keeps in the notation.
Therefore, when comparing matrix RPA susceptibilites to TPRF ones this flip of indices has to be
taken into account.

Mapping between spin-dependent and independent quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While the greek indices only carry orbital information, the latin indices used for the
quantities in TPRF carry orbital and spin information.
To map the spin-dependent to independent quantities and vice versa, we introduce the following notation

.. math::
    a = \mathrm{orb}(a)_{\sigma(a)} = \alpha_{\sigma(a)}\,,

where the :math:`\mathrm{orb}` function extracts the orbital information, mapping to the greek letters,
and the :math:`\sigma` extracts the spin.

With this we can state the mapping between the susceptibilites as

.. math::
    \chi^{(s)}(\bar{\alpha}\beta\bar{\gamma}\delta) =
    \chi^{\mathrm{RPA}}(\bar{\alpha}_\uparrow \beta_\uparrow \bar{\gamma}_\uparrow \delta_\uparrow)-
    \chi^{\mathrm{RPA}}(\bar{\alpha}_\uparrow \beta_\uparrow \bar{\gamma}_\downarrow \delta_\downarrow)\,,

.. math::
    \chi^{(c)}(\bar{\alpha}\beta\bar{\gamma}\delta) =
    \chi^{\mathrm{RPA}}(\bar{\alpha}_\uparrow \beta_\uparrow \bar{\gamma}_\uparrow \delta_\uparrow)+
    \chi^{\mathrm{RPA}}(\bar{\alpha}_\uparrow \beta_\uparrow \bar{\gamma}_\downarrow \delta_\downarrow)\,,

and

.. math::
    \chi^{\mathrm{RPA}}(\bar{a}b\bar{c}d) =
    \begin{cases}
    \frac{1}{2}\big(\chi^{(\mathrm{c})} + \chi^{(\mathrm{s})}\big)(\bar{\alpha}\beta\bar{\gamma}\delta),&
    \mathrm{if}\; \sigma(\bar{a}) = \sigma(b) = \sigma(\bar{c}) = \sigma(d)\\
    \frac{1}{2}\big(\chi^{(\mathrm{c})} - \chi^{(\mathrm{s})}\big)(\bar{\alpha}\beta\bar{\gamma}\delta),&
    \mathrm{if}\; \sigma(\bar{a}) = \sigma(b) \neq \sigma(\bar{c}) = \sigma(d)\\
    \chi^{(\mathrm{s})}(\bar{\alpha}\beta\bar{\gamma}\delta),&
    \mathrm{if}\; \sigma(\bar{a}) = \sigma(d) \neq \sigma(b) = \sigma(\bar{c})\\
    0, & \mathrm{else}
    \end{cases}\,.

And for the interaction they are given by

.. math::
    U^{(s)}(\alpha\bar{\beta}\gamma\bar{\delta}) =
    U(\alpha_\uparrow \bar{\beta}_\uparrow \gamma_\uparrow \bar{\delta}_\uparrow)-
    U(\alpha_\uparrow \bar{\beta}_\uparrow \gamma_\downarrow \bar{\delta}_\downarrow)\,,

.. math::
    U^{(c)}(\alpha\bar{\beta}\gamma\bar{\delta}) =
    -U(\alpha_\uparrow \bar{\beta}_\uparrow \gamma_\uparrow \bar{\delta}_\uparrow)-
    U(\alpha_\uparrow \bar{\beta}_\uparrow \gamma_\downarrow \bar{\delta}_\downarrow)\,,

.. math::
    U(a\bar{b}c\bar{d}) =
    \begin{cases}
    \frac{1}{2}\big(U^{(\mathrm{s})} - U^{(\mathrm{c})}\big)(\alpha\bar{\beta}\gamma\bar{\delta}),&
    \mathrm{if}\; \sigma(a) = \sigma(\bar{b}) = \sigma(c) = \sigma(\bar{d})\\
    \frac{1}{2}\big(-U^{(\mathrm{c})} - U^{(\mathrm{s})}\big)(\alpha\bar{\beta}\gamma\bar{\delta}),&
    \mathrm{if}\; \sigma(a) = \sigma(\bar{b}) \neq \sigma(c) = \sigma(\bar{d})\\
    U^{(\mathrm{s})}(\alpha\bar{\beta}\gamma\bar{\delta}),&
    \mathrm{if}\; \sigma(a) = \sigma(\bar{d}) \neq \sigma(\bar{b}) = \sigma(c)\\
    0, & \mathrm{else}
    \end{cases}\,.

Example
^^^^^^^

If you have a spin-dependent bare particle-hole bubble :samp:`chi00_wk` and a spin-dependent vertex
:samp:`U_abcd`, you could use the following code snippet to produce the corresponding
spin-dependent general RPA susceptibility :samp:`chi_wk`, without doing a spin-dependent calculation.

.. code-block:: python

 from triqs_tprf.rpa_tensor import lose_spin_degree_of_freedom
 chi00_wk_wo_spin = lose_spin_degree_of_freedom(chi00_wk, spin_fast=False)

 from triqs_tprf.rpa_tensor import lose_spin_degree_of_freedom
 U_c, U_s = split_quartic_tensor_in_charge_and_spin(U_abcd)

 from triqs_tprf.lattice import solve_rpa_PH
 chi_s = solve_rpa_PH(chi00_wk_wo_spin, U_s)
 chi_c = solve_rpa_PH(chi00_wk_wo_spin, -U_c) # Minus for correct charge rpa equation

 from triqs_tprf.rpa_tensor import general_susceptibility_from_charge_and_spin
 chi_wk = general_susceptibility_from_charge_and_spin(chi_c, chi_s, spin_fast=False)

Or you could already start at the spin-dependent Green's function :samp:`g0_wk` to construct
a spin-independent bare particle-hole bubble.

.. code-block:: python
 
 from triqs_tprf.rpa_tensor import lose_spin_degree_of_freedom
 g0_wk_wo_spin = lose_spin_degree_of_freedom(g0_wk, spin_fast=False)

 from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
 chi00_wk_wo_spin = imtime_bubble_chi0_wk(g0_wk_wo_spin, nw=1)


