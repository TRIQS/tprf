.. _eliashberg:

Linearized Eliashberg Equation
==============================

The linearized Eliashberg equation is a generalization of the linearized
Bardeen-Cooper-Schrieffer (BCS) gap equation to frequency dependent gaps.
It can be used to determine the critical (inverse) temperature
:math:`T_\mathrm{c}/\beta_\mathrm{c}`,
at which a transition to a superconducting state occurs,
and the symmetry of the corresponding superconducting gap function
:math:`\Delta^{\mathrm{s/t}}`.
It is given by

.. math::
    \Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta_\mathrm{c}}\sum_{K'}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(Q=0, K', K)
    G_{c\bar{f}}(K')G_{d\bar{e}}(-K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,.
    :label: linearized_eliashberg_1

where :math:`Q/K` is a combination of bosonic/fermionic Matsubara :math:`i\omega_n/i\nu_n` 
frequency and momentum :math:`\mathbf{k}`,
:math:`N_{\mathbf{k}}` is the number of momentum points,
:math:`\Gamma^{\mathrm{s/t}}` is the irreducible particle-particle vertex
and :math:`G` is the one-particle Green's function.

.. note::
   The bosonic Matsubara frequency and momentum in the particle-particle vertex is set to zero.
   This is because we are interested in Cooper-pairs which have a zero
   transfered momentum-frequency in a scattering process [#nourafkan]_.

.. note::
    The current implementation is restricted to :math:`SU(2)` symmetric systems.
    All indices are purely orbital and superconducting gaps :math:`\Delta` and
    particle-particle vertices :math:`\Gamma` are restricted to the singlet/triplet
    channel, shown by the superscripts s/t respectively. 
    But note, that the equations still hold for the spin-dependent case and
    one would soley need to implement the spin-dependent particle-particle vertex
    to use them.

Deriving the linearized Eliashberg equation from the normal state
-----------------------------------------------------------------

The singlet and triplet susceptibilties are given by

.. math::
    \chi^{\mathrm{s}}
    =
    -
    \chi^{(0), \mathrm{PP}}
    +
    \frac{1}{2}
    \chi^{(0), \mathrm{PP}}
    \mathbf{\Gamma}^{\mathrm{s}}
    \left[
    -
    \chi^{\mathrm{s}}
    +
    \chi^{(0), \mathrm{PP}}
    \right]
    \,,

and

.. math::
    \chi^{\mathrm{t}}
    =
    \chi^{(0), \mathrm{PP}}
    +
    \frac{1}{2}
    \chi^{(0), \mathrm{PP}}
    \mathbf{\Gamma}^{\mathrm{t}}
    \left[
    \chi^{\mathrm{t}}
    +
    \chi^{(0), \mathrm{PP}}
    \right]
    \,.

A transition from the normal state to a singlet/triplet superconducting one occurs
when the susceptibilties diverge.
This is the case, when the largest eigenvalue of 
:math:`\mp \frac{1}{2}\mathbf{\Gamma^{\mathrm{s/t}}} \mathbf{\chi}^{(0),{PP}}` becomes unity.
For a largest eigenvalues that is smaller than :math:`1` we are still in the
normal state,
but we can calculate the corresponding eigenvectors :math:`\Delta^{\mathrm{s/t}}`.
This corresponds to the following eigenvalue equation

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)
    = 
    \frac{1}{2N_{\mathbf{k}}^2 \beta^2}\sum_{K', K''}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(Q=0, K', K)
    \chi^{(0),{PP}}_{\bar{f}d\bar{e}c}(Q=0, K'', K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K'')\,,
 
where we incorporate the minus sign of the singlet channel in our definition of the 
singlet irreducible vertex to only keep track of one version of the Eliashberg equation.

We can write this like Eq. :eq:`linearized_eliashberg_1` by using the definiton
of :math:`\chi^{(0),{PP}}`

.. math::
    \chi^{(0),{PP}}_{\bar{a}b\bar{c}d}(Q=0, K, K') 
    =
    -N_{\mathbf{k}} \beta
    G_{d\bar{a}}(K)G_{b\bar{c}}(-K')\delta_{K, K'}\,,
    :label: chi_0_pp

which yields

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(Q=0, K', K)
    G_{c\bar{f}}(K')G_{d\bar{e}}(-K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,.
    :label: linearized_eliashberg_3

.. note::
    Our definiton of :math:`\chi^{(0),{PP}}` is different from [#bickers]_
    and [#nourafkan]_. This stems from the fact, that due to the indistinguishability 
    of the particles in the particle-particle channel doublecounting diagrams in the 
    Bethe-Salpeter equation (BSE) must be avoided. 
    We do this by defining the particle-particle BSE with a factor of
    :math:`\frac{1}{2}`, see :ref:`vertex` Eq. :eq:`BSE_PP`.
    In [#bickers]_ and [#nourafkan]_ the particle-particle BSE is defined without this
    factor and they include it in their definiton of :math:`\chi^{(0),{PP}}`.

This equation is valid for :math:`\lambda \leq 1`
and yields eigenvectors, which correspond to superconducting gap functions
that have not manifested yet.
At :math:`\lambda=1` the normal state breaks down and the superconducting
state with the corresponding gap emerges.
The size of the eigenvalues is therefore an indicator of how likely the associated gap
is to manifest.

Relation to the BCS gap equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In BCS theory the particle-particle vertex is considered to be
constant in a specific frequency range, which corresponds to gaps with
the same dependence.
For this case the summation over fermionic Matsubara frequencies in the linearized
Eliashberg equation Eq. :eq:`linearized_eliashberg_1` can be done analytically.
For a one-band case and a non-interacting Green's function with dispersion relation
:math:`\epsilon`, this yields

.. math::
    \Delta^{\mathrm{s/t}}(\mathbf{k}) =  -\frac{1}{2 N_{\mathbf{k}}}\sum_{\mathbf{k'}}
    \Gamma^{\mathrm{s/t}}(\mathbf{q}=\mathbf{0}, \mathbf{k}', \mathbf{k})
    \frac{\tan(\epsilon(\mathbf{k'})\beta/2)}{2\epsilon(\mathbf{k'})}
    \Delta^{\mathrm{s/t}}(\mathbf{k'})\,,
    :label: linearized_eliashberg_4

which corresponds to the linearized BCS gap equation.
The non-linear BCS gap equation can be obtained from Eq. :eq:`linearized_eliashberg_4` 
by substituting :math:`\epsilon` with
:math:`\sqrt{\epsilon(\mathbf{k})^2 + |\Delta(\mathbf{k})|^2}`.


Details for applications 
------------------------

SPOT Condition
^^^^^^^^^^^^^^

In the general case the superconducting gap function :math:`\Delta` is dependent on 
momentum :math:`\mathbf{k}`, fermionic Matsubara frequency :math:`i\nu_n`,
orbital-indices :math:`a,b` and spin-indices :math:`\alpha,\beta`

.. math::
    \Delta \equiv \Delta_{a\alpha;b\beta}(i\nu, \mathbf{k})\,.

Because the Pauli principle dictates :math:`\Delta` to be odd under particle exchange,
the symmetry combinations of those four degrees of freedom are constrained.
This is formalized as the so called :math:`SPOT` condition

.. math::
    \hat{S}\hat{P}\hat{O}\hat{T} \Delta_{a\alpha;b\beta}(i\nu, \mathbf{k}) 
    =
    - \Delta_{b\beta;a\alpha}(-i\nu, -\mathbf{k})\,,
    
with the operators :math:`\hat{S}`, :math:`\hat{P}`, :math:`\hat{O}`, :math:`\hat{T}`,
that denote permutation of electrons in spin space (:math:`\hat{S}`),
real space (parity) (:math:`\hat{P}`),
orbital space (:math:`\hat{O}`), and time (frequency) (:math:`\hat{T}`).
While :math:`\Delta` has to be odd under the combined action of the symmetry operations
:math:`\hat{S}\hat{P}\hat{O}\hat{T}`,
it can either be even (:math:`+`) or odd (:math:`-`) under each separate operation,
i.e.

.. math::
    \hat{S}\Delta_{a\alpha;b\beta}(i\nu, \mathbf{k}) 
	&=
	\pm \Delta_{a\beta;b\alpha}(i\nu, \mathbf{k})\,,\\
	\hat{P}\Delta_{a\alpha;b\beta}(i\nu, \mathbf{k}) 
	&=
	\pm \Delta_{a\alpha;b\beta}(i\nu, -\mathbf{k})\,,\\
	\hat{O}\Delta_{a\alpha;b\beta}(i\nu, \mathbf{k}) 
	&=
	\pm \Delta_{b\alpha;a\beta}(i\nu, \mathbf{k})\,,\\
	\hat{T}\Delta_{a\alpha;b\beta}(i\nu, \mathbf{k}) 
	&=
	\pm \Delta_{a\alpha;b\beta}(-i\nu, \mathbf{k})\,.

A gap function can therefore be classified as even (:math:`+`) or odd (:math:`-`)
under these four degrees of freedom. By calculating the superconducting gap in the
singlet/triplet channel, we fix the spin symmetry to odd/even respectively.
This leaves us with four symmetry combinations for both singlet and triplet gaps,
which we list in the table below.

.. table:: 
    :align: center

    +-----------------------------------------------+-----------------------------------------------+
    |                  Spin-singlet                 |                  Spin-triplet                 |
    +===========+===========+===========+===========+===========+===========+===========+===========+
    |     S     |     P     |     O     |     T     |     S     |     P     |     O     |     T     |
    +-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
    | :math:`-` | :math:`+` | :math:`+` | :math:`+` | :math:`+` | :math:`-` | :math:`-` | :math:`-` |
    +-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
    | :math:`-` | :math:`-` | :math:`-` | :math:`+` | :math:`+` | :math:`+` | :math:`+` | :math:`-` |
    +-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
    | :math:`-` | :math:`-` | :math:`+` | :math:`-` | :math:`+` | :math:`+` | :math:`-` | :math:`+` |
    +-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+
    | :math:`-` | :math:`+` | :math:`-` | :math:`-` | :math:`+` | :math:`-` | :math:`+` | :math:`+` |
    +-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+

Because all other combinations are unphysical it is possible to restrict the gap to the
allowed symmetries while solving the linearized Eliashberg equation. 

.. _eliashberg_rpa:

Random phase approximation for the irreducible particle-particle vertex
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The irreducible particle-particle vertex is given by the parquet equation,
which can be expressed in terms of the fully irreducible vertex :math:`\Lambda`
and the channel reducible vertex ladder functions :math:`\Phi`.
It is given in the singlet channel by

.. math::
    \Gamma^{\text{s}}_{a\overline{b}c\overline{d}}(Q, K, K') =&
	-
	\Lambda^{\text{s}}_{a\overline{b}c\overline{d}}(Q, K, K')
	+
	\left[
	\frac{3}{2}
	\Phi^{\text{m}}_{a\overline{b}c\overline{d}}
	-
	\frac{1}{2}
	\Phi^{\text{d}}_{a\overline{b}c\overline{d}}
	\right](Q-K-K', K, K')
    \\
	&+
	\left[
	\frac{3}{2}
	\Phi^{\text{m}}_{c\overline{b}a\overline{d}}
	-
	\frac{1}{2}
	\Phi^{\text{d}}_{c\overline{b}a\overline{d}}
	\right](K-K', Q-K, K')
    :label: singlet_gamma_no_approx

and in the triplet channel by

.. math::
    \Gamma^{\text{t}}_{a\overline{b}c\overline{d}}(Q, K, K') =&
    \Lambda^{\text{t}}_{a\overline{b}c\overline{d}}(Q, K, K')
    +
    \left[
    \frac{1}{2}
    \Phi^{\text{m}}_{a\overline{b}c\overline{d}}
    +
    \frac{1}{2}
    \Phi^{\text{d}}_{a\overline{b}c\overline{d}}
    \right](Q-K-K', K, K')
    \\
    &+
    \left[
    -
    \frac{1}{2}
    \Phi^{\text{m}}_{c\overline{b}a\overline{d}}
    -
    \frac{1}{2}
    \Phi^{\text{d}}_{c\overline{b}a\overline{d}}
    \right](K-K', Q-K, K')
    \,,
    :label: triplet_gamma_no_approx

with the spin diagonalized reducible vertex ladder functions given by

.. math::
    \Phi^{\text{d/m}}_{a\overline{b}c\overline{d}}(Q, K, K')
    =
    \frac{1}{(N_\mathbf{k}\beta)^2}
    \sum_{K'', K'''}
    \Gamma^{\text{d/m}}(Q, K, K'') \chi^{\text{d/m}}(Q, K'', K''') \Gamma^{\text{d/m}}(Q, K''', K')
    \,.

Note, that the superscripts :math:`\mathrm{d/m}` indicate the density/magnetic channel.

Now, in the random phase approximation (RPA) the susceptibilities :math:`\chi^{\text{d/m}}`
are approximated by the RPA bubble susceptibility,
and all vertices are substituted by the local and static bare Kanamori interaction :math:`U^{\mathrm{d/m}}`,
given by

.. math::
    U^{\mathrm{d/m}}_{a\bar{b}c\bar{d}} =
    \begin{cases}
    U/U, & \mathrm{if}\;a=\bar{b}=c=\bar{d} \\
    -U'+2J/U', & \mathrm{if}\;a=\bar{d}\neq \bar{b}=c \\
    2U'-J/J, & \mathrm{if}\;a=\bar{b}\neq c=\bar{d} \\
    J/J, & \mathrm{if}\;a=c\neq \bar{b}=\bar{d} \\
    0, & \mathrm{else}
    \end{cases}\,,

with the Hubbard interaction :math:`U` and the Hund's :math:`J`.
The reducible ladder vertices then beceome only dependent on one bosonic Frequence and
momentum pair :math:`Q`

.. math::
    \Phi^{\text{d/m}}_{a\overline{b}c\overline{d}}(Q)
    &\approx
    \frac{1}{(N_\mathbf{k}\beta)^2}
    \sum_{K'', K'''}
    U^{\text{d/m}}\chi^{\text{d/m}}(Q, K'', K''') U^{\text{d/m}}
    \\
    &\approx
    U^{\mathrm{d/m}}
    \chi^{\text{d/m}}(Q) U^{\mathrm{d/m}}
    \,,

and the fully irreducible vertices become

.. math::
    \Lambda^{\mathrm{s}}
    \approx
    -
    \frac{1}{2}U^{\mathrm{d}}
    -
    \frac{3}{2}U^{\mathrm{m}}
    \,,

.. math::
    \Lambda^{\mathrm{t}}
    \approx
    -
    \frac{1}{2}U^{\mathrm{d}}
    +
    \frac{1}{2}U^{\mathrm{m}}
    \,.

In this approximation the irreducible singlet/triplet vertex for :math:`Q=0` takes the form

.. math::
    \Gamma^{\text{s}}_{a\overline{b}c\overline{d}}(Q=0, K, K') =&
	\frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{d}}
	+
	\frac{3}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{m}}
	+
	\left[
	\frac{3}{2}
	\Phi^{\text{m}}_{a\overline{b}c\overline{d}}
	-
	\frac{1}{2}
	\Phi^{\text{d}}_{a\overline{b}c\overline{d}}
	\right](-K-K')
    \\
	&+
	\left[
	\frac{3}{2}
	\Phi^{\text{m}}_{c\overline{b}a\overline{d}}
	-
	\frac{1}{2}
	\Phi^{\text{d}}_{c\overline{b}a\overline{d}}
	\right](K-K')
	\,,
    :label: singlet_gamma

and

.. math::
    \Gamma^{\text{t}}_{a\overline{b}c\overline{d}}(Q=0, K, K') =&
	-
	\frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{d}}
	+
	\frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{m}}
	+
	\left[
	\frac{1}{2}
	\Phi^{\text{m}}_{a\overline{b}c\overline{d}}
	+
	\frac{1}{2}
	\Phi^{\text{d}}_{a\overline{b}c\overline{d}}
	\right](-K-K')
    \\
	&+
	\left[
	-
	\frac{1}{2}
	\Phi^{\text{m}}_{c\overline{b}a\overline{d}}
	-
	\frac{1}{2}
	\Phi^{\text{d}}_{c\overline{b}a\overline{d}}
	\right](K-K')
    \,.
    :label: triplet_gamma

Note, that inserting both the singlet :eq:`singlet_gamma` and the triplet vertex 
:eq:`triplet_gamma` in the linearized Eliashberg equation :eq:`linearized_eliashberg_1`
flips the momentum/frequency dependence. We therefore get inside the Eq. :eq:`linearized_eliashberg_1` 
density and magnetic vertices :math:`\Phi^{\text{d/m}}` with an index flip and with a :math:`K'-K`
dependence, :math:`\Phi_{c\overline{b}a\overline{d}}(K'-K)`, and without an index flip 
and a :math:`-K'-K` dependence, :math:`\Phi_{a\overline{b}c\overline{d}}(-K'-K)`.
These two terms can be transformed into each other by abiding the frequency, momentum and orbital
symmetry of the gap. 
For example :math:`\Phi_{a\overline{b}c\overline{d}}(-K'-K)` transforms into 
:math:`\pm\Phi_{c\overline{b}a\overline{d}}(K'-K)` for a singlet/triplet gap.
By using the property

.. math::
   \Phi^{\text{d/m}}(K'-K)
   =
   \text{Complex Conjugate}(\Phi^{\text{d/m}}(K-K'))

for a translation symmetric system we can write Eq. :eq:`singlet_gamma` and :eq:`triplet_gamma` as

.. math::
    \Gamma^{\text{s}}_{a\overline{b}c\overline{d}}(K-K') \equiv
	\frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{d}}
	+
	\frac{3}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{m}}
	+
	\text{Complex Conjugate}
	\left[
	3 
	\Phi^{\text{m}}_{c\overline{b}a\overline{d}}(K-K')
	-
	\Phi^{\text{d}}_{c\overline{b}a\overline{d}}(K-K')
	\right]
	\,,
    :label: singlet_gamma_2

.. math::
    \Gamma^{\text{t}}_{a\overline{b}c\overline{d}}(K - K') \equiv
	-
	\frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{d}}
	+
	\frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{m}} 
	+
	\text{Complex Conjugate}
	\left[
    -
	\Phi^{\text{m}}_{c\overline{b}a\overline{d}}(K-K')
	-
	\Phi^{\text{d}}_{c\overline{b}a\overline{d}}(K-K')
	\right]
	\,.
    :label: triplet_gamma_2

Note, that this simplification is only allowed if the solutions of :math:`\Delta^{\mathrm{s/t}}`
are restricted to the allowed symmetries, otherwise unphysical solutions can occur.
Also note, that the RPA particle-particle vertices in
Eq. :eq:`singlet_gamma_2` and :eq:`triplet_gamma_2` only depend on the difference
between the two fermionic Matsubara frequencies, i.e. a bosonic Matsubara frequency and one momentum.
We can therefore write the linearized Eliashberg equation
:eq:`linearized_eliashberg_3` as

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(K-K')
    G_{c\bar{f}}(K')G_{d\bar{e}}(-K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,,
    :label: linearized_eliashberg_5

which is the **form it is implemented as now** in :meth:`triqs_tprf.eliashberg.solve_eliashberg`.
    
This allows us to get rid of the summation by using the convolution theorem

.. math::
    \lambda
    \mathcal{F}\left[\Delta_{\bar{a}\bar{b}}^{\mathrm{s/t}}(K)\right]=  -\frac{1}{2}
    \mathcal{F}\left[\Gamma_{c\bar{a}d\bar{b}}^{\mathrm{s/t}}(K-K')\right]
    \mathcal{F}\left[
    G_{c\bar{f}}(K')G_{d\bar{e}}(-K')
    \Delta_{\bar{e}\bar{f}}^{\mathrm{s/t}}(K')
    \right]\,,
    :label: linearized_eliashberg_5

making the calculation computationaly more efficient for large numbers of frequencies
and momenta.
But note, that for small numbers of frequencies and/or momenta using the sum
instead of the convolution theorem can be more effecient.

.. note::
    It is possible to expand the current implementation of the Eliashberg equation to
    also allow for irreducible vertices to be explicitly dependent on two fermionic
    frequency and momenta pairs.
    For an idea on how to tackle such a task see the following draft
    `here <https://github.com/TRIQS/tprf/blob/eliashberg_with_phi/c%2B%2B/triqs_tprf/lattice/eliashberg.cpp#L289>`_
    and 
    `here <https://github.com/TRIQS/tprf/blob/eliashberg_with_phi/python/triqs_tprf/eliashberg.py#L216>`_.


.. rubric:: References

.. [#abrikosov] A. A. Abrikosov, L. P. Gor’kov, and I. E. Dzyaloshinski, Pergamon, Oxford (1965)
.. [#yanase] Y. Yanase, T. Jujo, T. Nomura, et. al., Physics Reports 387, 1-149 (2003)
.. [#takimoto] T. Takimoto, T. Hotta, and K. Ueda, PRB 69, 104504 (2004)
.. [#bickers] N. E. Bickers, Self-Consistent Many-Body Theory for Condensed Matter Systems. Theoretical Methods for Strongly Correlated Electrons, 237–296. 6 (2006)
.. [#rohringer] G. Rohringer, New routes towards a theoretical treatment of nonlocal electronic correlations (2013)
.. [#nourafkan] R. Nourafkan, G. Kotliar, and A. M. Tremblay, Physical Review Letters 117, 1, (Supplementary) (2016) 
.. [#linder] J. Linder and A. V. Balatsky, Reviews of Modern Physics 91, 45005 (2019)
