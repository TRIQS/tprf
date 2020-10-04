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
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
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
   transfered momentum-frequency in a scattering process.

.. note::
    The current implementation is restricted to :math:`SU(2)` symmetric systems.
    All indices are purely orbital and superconducting gaps :math:`\Delta` and
    particle-particle vertices :math:`\Gamma` are restricted to the singlet/triplet
    channel, shown by the superscripts s/t respectively. 

Deriving the linearized Eliashberg equation from the normal state
-----------------------------------------------------------------

Generally speaking a transition from the normal state to the superconducting
one occurs when the particle-particle susceptibility diverges.

.. math::
    \mathbf{\chi}^{\mathrm{s/t}} = [\mathbf{1}-\mathbf{\Gamma}^{\mathrm{s/t}}
    \mathbf{\chi}^{(0),{PP}}]^{-1}
    \mathbf{\chi}^{(0),{PP}}

This is the case when the largest eigenvalue of 
:math:`\mathbf{\Gamma^{\mathrm{s/t}}} \mathbf{\chi}^{(0),{PP}}` becomes unity.
For a largest eigenvalues that is smaller than :math:`1` we are still in the
normal state,
but we can calculate the corresponding eigenvectors :math:`\Delta^{\mathrm{s/t}}`.
This corresponds to the following eigenvalue equation

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)
    = 
    \frac{1}{N_{\mathbf{k}}^2 \beta^2}\sum_{K', K''}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    \chi^{(0),{PP}}_{\bar{e}d\bar{f}c}(Q=0, K', K'')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,,
    :label: linearized_eliashberg_2

which we can write like Eq. :eq:`linearized_eliashberg_1` with the definiton
of :math:`\chi^{(0),{PP}}`

.. math::
    \chi^{(0),{PP}}_{\bar{a}b\bar{c}d}(Q, K, K') 
    =
    -\frac{N_{\mathbf{k}} \beta}{2}
    G_{d\bar{a}}(K)G_{b\bar{c}}(-K')\delta_{K, K'}\,,
    :label: chi_0_pp

as

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,.
    :label: linearized_eliashberg_3

.. note::
    There is an inconsistency with a factor of :math:`\frac{1}{2}` with
    Eq. :eq:`chi_0_pp` and Eq. :eq:`bare_pp_sus_def`.
    As there is no bare particle-particle bubble implementation yet,
    and the Eliashberg implementation is self-consistent,
    we don't have any problems. 
    But for future implementations this needs to be addressed.

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
    \Gamma^{\mathrm{s/t}}(\mathbf{q}=\mathbf{0}, \mathbf{k}, \mathbf{k'})
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
and the channel reducible vertex-ladder functions :math:`\Phi`.
It is given in the singlet channel by

.. math::
    \Gamma^{\mathrm{s}}_{a\bar{b}c\bar{d}}(Q=0, K, K')
	\equiv &
	\frac{3}{2}
	\left[
	\Phi^{\mathrm{m}}_{a\bar{b}c\bar{d}}(K-K')
	+
	\Phi^{\mathrm{m}}_{c\bar{b}a\bar{d}}(K+K')
	\right]
	\\&-
	\frac{1}{2}
	\left[
	\Phi^{\mathrm{d}}_{a\bar{b}c\bar{d}}(K-K')
	+
	\Phi^{\mathrm{d}}_{c\bar{b}a\bar{d}}(K+K')
	\right]
	+
	\Lambda^{\mathrm{s}}_{a\bar{b}c\bar{d}}\,,
    :label: singlet_gamma

and in the triplet channel by

.. math::
    \Gamma^{\mathrm{t}}_{a\bar{b}c\bar{d}}(Q=0, K, K')
	\equiv &
	-\frac{1}{2}
	\left[
	\Phi^{\mathrm{m}}_{a\bar{b}c\bar{d}}(K-K')
	-
	\Phi^{\mathrm{m}}_{c\bar{b}a\bar{d}}(K+K')
	\right]
	\\&-
	\frac{1}{2}
	\left[
	\Phi^{\mathrm{d}}_{a\bar{b}c\bar{d}}(K-K')
	-
	\Phi^{\mathrm{d}}_{c\bar{b}a\bar{d}}(K+K')
	\right]
	+
	\Lambda^{\mathrm{t}}_{a\bar{b}c\bar{d}}\,,
    :label: triplet_gamma


where the vertex-ladder functions are given by

.. math::
    \Phi^{\text{d/m}}_{a\overline{b}c\overline{d}}(Q)
    =
    \Lambda^{\text{d/m}} \chi^{\text{d/m}}(Q) \Lambda^{\text{d/m}}\,.


Note, that the superscripts :math:`\mathrm{d/m}` indicate the density/magnetic channel.

Now, in the random phase approximation (RPA) the susceptibilities :math:`\chi^{\text{d/m}}`
are approximated by the RPA bubble susceptibility,
and the vertices are approximated by

.. math::
    \Lambda^{\text{d/m}} \approx U^{\mathrm{d/m}}\,,

and

.. math::
    \Lambda^{\text{s/t}} \approx \frac{1}{2}(U^{\mathrm{d}} + U^{\mathrm{m}})\,.

Here :math:`U^{\mathrm{d/m}}` is the bare local Kanamori interaction given by 

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

Note, that in both singlet :eq:`singlet_gamma` and
triplet :eq:`triplet_gamma` a density and magnetic
:math:`\Phi` term appears twice.
Once without an index flip and a dependence on :math:`K-K'`,
:math:`\Phi_{a\overline{b}c\overline{d}}(K-K')`,
and another time with an index flip and a dependence on :math:`K+K'`, 
:math:`\Phi_{c\overline{b}a\overline{d}}(K+K')`.
Inside the linearized Eliashberg equation :eq:`linearized_eliashberg_3`
the :math:`\Phi_{c\overline{b}a\overline{d}}(K+K')` term
picks up a sign which depends on the frequency, momentum and orbital
symmetry of the gap :math:`\Delta^{\mathrm{s/t}}`.
For all allowed singlet combinations it is positive and for all allowed triplet ones
negative. Therefore Eq. :eq:`singlet_gamma` and Eq. :eq:`triplet_gamma` become

.. math::
    \Gamma^{\text{s}}_{a\overline{b}c\overline{d}}(Q=0, K, K') \equiv
    3 
    \Phi^{\text{m}}_{a\overline{b}c\overline{d}}(K-K')
    -
    \Phi^{\text{d}}_{a\overline{b}c\overline{d}}(K-K')
    +
    \Lambda^{\text{s}}_{a\overline{b}c\overline{d}}
    \,,
    :label: singlet_gamma_2

.. math::
    \Gamma^{\text{t}}_{a\overline{b}c\overline{d}}(Q=0, K, K') \equiv
    -
    \Phi^{\text{m}}_{a\overline{b}c\overline{d}}(K-K')
    -
    \Phi^{\text{d}}_{a\overline{b}c\overline{d}}(K-K')
    +
    \Lambda^{\text{t}}_{a\overline{b}c\overline{d}} 
    \,.
    :label: triplet_gamma_2

Note, that this simplification is only allowed, if the solutions of :math:`\Delta^{\mathrm{s/t}}`
are restricted to the allowed symmetries, otherwise unphysical solution can occur.
Also note, that the RPA particle-particle vertices in
Eq. :eq:`singlet_gamma_2` and :eq:`triplet_gamma_2` only depend on the difference
between the two fermionic Matsubara frequencies, i.e. a bosonic Matsubara frequency and one momentum.
We can therefore write the linearized Eliashberg equation
:eq:`linearized_eliashberg_3` as

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(K-K')
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,,
    :label: linearized_eliashberg_5

which is the **form it is implemented as now** in :meth:`triqs_tprf.eliashberg.solve_eliashberg`.
    
This allows us to get rid of the summation by using the convolution theorem

.. math::
    \lambda
    \mathcal{F}\left[\Delta_{\bar{a}\bar{b}}^{\mathrm{s/t}}(K)\right]=  -\frac{1}{2}
    \mathcal{F}\left[\Gamma_{c\bar{a}d\bar{b}}^{\mathrm{s/t}}(K-K')\right]
    \mathcal{F}\left[
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
    \Delta_{\bar{e}\bar{f}}^{\mathrm{s/t}}(K')
    \right]\,,
    :label: linearized_eliashberg_5

making the calculation computationaly more efficient.

.. rubric:: References

.. [#abrikosov] A. A. Abrikosov, L. P. Gor’kov, and I. E. Dzyaloshinski, Pergamon, Oxford (1965)
.. [#yanase] Y. Yanase, T. Jujo, T. Nomura, et. al., Physics Reports 387, 1-149 (2003)
.. [#takimoto] T. Takimoto, T. Hotta, and K. Ueda, PRB 69, 104504 (2004)
.. [#bickers] N. E. Bickers, Self-Consistent Many-Body Theory for Condensed Matter Systems. Theoretical Methods for Strongly Correlated Electrons, 237–296. 6 (2006)
.. [#rohringer] G. Rohringer, New routes towards a theoretical treatment of nonlocal electronic correlations (2013)
.. [#nourafkan] R. Nourafkan, G. Kotliar, and A. M. Tremblay, Physical Review Letters 117, 1, (Supplementary) (2016) 
.. [#linder] J. Linder and A. V. Balatsky, Reviews of Modern Physics 91, 45005 (2019)
