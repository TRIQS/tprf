.. _eliashberg:

Linearized Eliashberg Equation
================================

The linearized Eliashberg equation is a generalization of the
Bardeen-Cooper-Schrieffer (BCS) gap equation to frequency dependent gaps.
It can be used to determine the critical temperature :math:`T_\mathrm{c}`,
at which a transition to a superconducting state occurs,
and the symmetry of the corresponding gap function :math:`\Delta`.
It is given by

.. math::
    \Delta(K) =  \frac{T^2_{\mathrm{c}}}{2 N_{\mathbf{k}}^2}\sum_{K', K''}
    \Gamma^{PP}(Q=0, K, K')
    \chi^{(0),{PP}}(Q=0, K', K'')
    \Delta(K'')\,,
    :label: linearized_eliashberg_1

where :math:`Q/K` is a combination of bosonic/fermionic Matsubara 
frequency and momentum,
:math:`N_{\mathbf{k}}` is the number of momentum points,
:math:`\Gamma^{PP}` is the irreducible particle-particle vertex,
:math:`\chi^{(0),{PP}}` is the bare particle-particle susceptibility,
as defined in :eq:`bare_pp_sus_def`,
and :math:`\Delta` is the gap.

Note, that the bosonic Matsubara frequency and momentum is set to zero.
This is because we are interested in Cooper-pairs which have a zero
transfered momentum-frequency in a scattering process.
    
Normal state and superconducting state
--------------------------------------

The linearized Eliashberg equation can be seen from two perspectives.
On one hand we are in the normal state and want to find the transition
to the superconducting one,
and on the other we are in the superconducting state using the limit of
small gaps :math:`\Delta \ll 1`.

Normal state
^^^^^^^^^^^^

Generally speaking a transition from the normal state to a superconducting
one occurs when the particle-particle susceptibility diverges.
This is the case when the product of irreducible vertex and
bare susceptibility becomes unity

.. math::
    \Gamma^{PP} \cdot \chi^{(0),{PP}} = 1\,.

For values smaller than :math:`1` we are still in the normal state,
but we can calculate the corresponding eigenvectors :math:`\Delta`.
This corresponds to an eigenvalue equation of Eq. :eq:`linearized_eliashberg_1`,
which we can write with the definiton of :math:`\chi^{(0),{PP}}` 
in :eq:`bare_pp_sus_def` as,

.. math::
    \lambda\Delta(K) =  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{PP}(Q=0, K, K')
    G(K')G(-K')
    \Delta(K')\,.
    :label: linearized_eliashberg_2

The eigenvalue :math:`\lambda` is therefore an indicator of how strong
a superconducting instability is.
These calculations are valid for :math:`\lambda \leq 1`
and yield gaps which have not actually manifested yet.
At :math:`\lambda=1` the normal state breaks down and the superconducting
state with the corresponding gap emerges.

Superconducting state
^^^^^^^^^^^^^^^^^^^^^

For a calculation in the superconducting state we must extend the basis
to accommodate Cooper-pairs.
In addition to the normal single-particle Green's function

.. math::
   G_{a\bar{b}}(\tau, \mathbf{r}) 
   \equiv 
   - \langle \mathcal{T} c_{a}(\tau, \mathbf{r})
   c^\dagger_{\bar{b}}(0, \mathbf{0}) \rangle \,,

and its backwards propagating counterpart

.. math::
   \overline{G}_{\bar{a}b}(\tau, \mathbf{r}) 
   \equiv 
   - \langle \mathcal{T} c^\dagger_{\bar{a}}(\tau, \mathbf{r})
   c_{b}(0, \mathbf{0}) \rangle\,,

we have to introduce the single-particle anomalous Green's functions
:math:`F` and :math:`\overline{F}`.
These are defined as

.. math::
    F_{ab}(\tau, \mathbf{r}) 
    \equiv
   \langle \mathcal{T} c_{a}(\tau, \mathbf{r}) 
   c_{b}(0, \mathbf{0}) \rangle
   \,,

and

.. math::
   \overline{F}_{\bar{a}\bar{b}}(\tau, \mathbf{r}) 
   \equiv
   \langle \mathcal{T} c^\dagger_{\bar{a}}(\tau, \mathbf{r}) 
   c^\dagger_{\bar{b}}(0, \mathbf{0}) \rangle\,.

Fourier transforming to Matsubara frequency and momentum space then gives that

.. math::
   \overline{G}_{\bar{a}b}(i\nu_n, \mathbf{k}) = 
   -G_{b\bar{a}}(-i\nu_n, -\mathbf{k})\,,
   :label: g_bar_to_g

and

.. math::
   \overline{F}_{\bar{a}\bar{b}}(i\nu_n, \mathbf{k})
   =
   [F_{ab}(-i\nu_n, \mathbf{k}) ]^{*}\,.
   :label: f_bar_to_f

All four Green's functions are coupled and given by

.. math::
   \left( \begin{array}{cc}
     G & F \\
     \overline{F} & \overline{G}\\
   \end{array} \right)
   =
   \left( \begin{array}{cc}
     \left(G^{(0)}\right)^{-1} - \Sigma & \Delta \\
     \overline{\Delta} & \left(\overline{G}^{(0)}\right)^{-1} - \overline{\Sigma} \\
   \end{array} \right)^{-1}\,,

where :math:`\Sigma` and :math:`\overline{\Sigma}` are the normal self-energies
and :math:`\Delta` and :math:`\overline{\Delta}` the anomalous ones,
which correspond to the gap function.

The anomalous self-energy is given by

.. math::
    \Delta(K)
    =
    -\frac{1}{N_{\mathbf{k}} \beta} \sum_{K'} \Gamma^{PP}(Q=0, K, K') F(K')\,.
    :label: anomalous_self_energy

with

.. math::
    F(K)
    =
    \frac{\Delta(K)}
    {\left(\left(G^{(0)}(K)\right)^{-1}-\Sigma(K)\right)
    \left(\left(G^{(0)}(K)\right)^{-1}-\Sigma(K)\right)
    -
    \overline{\Delta}(K) \Delta(K)}\,.
    :label: f_explicit

With this we could calculate the gap for the superconducting state
below :math:`T_\mathrm{c}`,
but note, that this is not that trivial,
because the self-energy :math:`\Sigma` is also coupled to :math:`\Delta`.
In the limit of :math:`\Delta \ll 1`,
i.e. linearizing Eq. :eq:`f_explicit`,
we also decouple the self-energy and
Eq. :eq:`anomalous_self_energy` becomes the the linearized Eliashberg equation.

Relation to the BCS gap equation
----------------------------

In BCS theory the particle-particle pairing vertex is considered to be
constant in a specific frequency range, which corresponds to gaps with
the same dependence.
For this case the summation over fermionic Matsubara frequencies in
Eq. :eq:`linearized_eliashberg_2` can be done analytically

.. math::
    \lambda\Delta(\mathbf{k}) =  -\frac{1}{2 N_{\mathbf{k}}}\sum_{\mathbf{k'}}
    \Gamma^{PP}(\mathbf{q}=\mathbf{0}, \mathbf{k}, \mathbf{k'})
    \frac{\tan(\epsilon(\mathbf{k'})\beta/2)}{2\epsilon(\mathbf{k'})}
    \Delta(\mathbf{k'})\,.
    :label: linearized_eliashberg_3

Here we assumed a non-interacting Gree'ns function with dispersion relation
:math:`\epsilon`.
Eq. :eq:`linearized_eliashberg_3` corresponds to the linearized BCS gap
equation.
The BCS gap equation can be obtained from Eq. :eq:`linearized_eliashberg_3` 
by substituting :math:`\epsilon` with
:math:`\sqrt{\epsilon(\mathbf{k})^2 + |\Delta(\mathbf{k})|^2}`.
This is equivalent to using Eq. :eq:`anomalous_self_energy`
with a non-linearized anomalous Green's function :eq:`f_explicit`.

Spin diagonalization
--------------------

All objects above were formulated with combined spin and orbital indices.
For :math:`SU(2)` symmetric systems we can drop the spin dependency
by diagonalizing everything in spin.
This diagonalization splits the superconducting gap :math:`\Delta`
in two channels,
singlet

.. math::
    \Delta^{\mathrm{s}}
    =
    \Delta_{\uparrow\downarrow}
    -
    \Delta_{\downarrow\uparrow}\,,

and triplet

.. math::
    \Delta^{\mathrm{t}}
    =
    \Delta_{\uparrow\downarrow}
    +
    \Delta_{\downarrow\uparrow}\,.

We can then express Eq. :eq:`linearized_eliashberg_2` in either of those two
channels.
Doing this for the singlet channel,
while suppressing frequency, momentum and orbital indices, yields

.. math::
    \lambda
    \Delta^{\mathrm{s}}
    &=
    \lambda
    \left(\Delta_{\uparrow\downarrow} - \Delta_{\downarrow\uparrow}\right)
    \\
    &=
    -\left[\Gamma^{PP}_{\uparrow\uparrow\downarrow\downarrow}
    \;
    G_{\uparrow\uparrow}G_{\downarrow\downarrow}
    \Delta_{\uparrow\downarrow}
    +
    \Gamma^{PP}_{\uparrow\downarrow\downarrow\uparrow}
    \;
    G_{\downarrow\downarrow}G_{\uparrow\uparrow}
    \Delta_{\downarrow\uparrow}
    \right]
    +
    \\
    &\quad\quad
    \left[
    \Gamma^{PP}_{\downarrow\downarrow\uparrow\uparrow}
    \;
    G_{\downarrow\downarrow}G_{\uparrow\uparrow}
    \Delta_{\downarrow\uparrow}
    -
    \Gamma^{PP}_{\downarrow\uparrow\uparrow\downarrow}
    \;
    G_{\uparrow\uparrow}
    G_{\downarrow\downarrow}
    \Delta_{\uparrow\downarrow}
    \right]
    \\
    &=
    -\Gamma^{PP}_{\uparrow\uparrow\downarrow\downarrow}
    \;
    GG
    \left(
    \Delta_{\uparrow\downarrow}
    -
    \Delta_{\downarrow\uparrow}
    \right)
    +
    \Gamma^{PP}_{\uparrow\downarrow\downarrow\uparrow}
    \;
    GG
    \left(
    \Delta_{\uparrow\downarrow}
    -
    \Delta_{\downarrow\uparrow}
    \right)
    \\
    &=
    -\left(\Gamma^{PP}_{\uparrow\uparrow\downarrow\downarrow}
    -\Gamma^{PP}_{\uparrow\downarrow\downarrow\uparrow}
    \;
    \right)
    GG
    \Delta^{\mathrm{s}}
    \\
    &=
    -\Gamma^{\mathrm{s}}
    GG
    \Delta^{\mathrm{s}}\,.

This is analog for the triplet channel and we obtain for the spin diagonalized
linearized Eliashberg equation

.. math::
    \lambda\Delta^{\mathrm{s/t}}(K) =  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}(Q=0, K, K')
    G(K')G(-K')
    \Delta^{\mathrm{s/t}}(K')\,.
    :label: linearized_eliashberg_4

Random phase approximation for the irreducible particle-particle vertex
-----------------------------------------------------------------------

To obtain the irreducible particle-particle vertex in the RPA
one substitutes all vertices with the bare one in the parquet equation.
In the spin diagonalized form this yields for the singlet

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
	U^{\mathrm{s}}_{a\bar{b}c\bar{d}}\,,
    :label: singlet_gamma

and for the triplet

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
	U^{\mathrm{s}}_{a\bar{b}c\bar{d}}\,,
    :label: triplet_gamma

with

.. math::
	\Phi^{\mathrm{d/m}}_{a\bar{b}c\bar{d}}(Q)
	\equiv
	U^{\mathrm{d/m}}\chi^{\mathrm{d/m}}(Q)U^{\mathrm{d/m}}\,.

Note, that the superscripts :math:`\mathrm{d}` and :math:`\mathrm{m}`
indicate the density and magnetic channel.
The bare vertex in its respective channel is given by

.. math::
    U^{\mathrm{d/m}}_{a\bar{b}c\bar{d}} =
    \begin{cases}
    U/U, & \mathrm{if}\;a=\bar{b}=c=\bar{d} \\
    -U'+2J/U', & \mathrm{if}\;a=\bar{d}\neq \bar{b}=c \\
    2U'-J/J, & \mathrm{if}\;a=\bar{b}\neq c=\bar{d} \\
    J/J, & \mathrm{if}\;a=c\neq \bar{b}=\bar{d} \\
    0, & \mathrm{else}
    \end{cases}\,.

For an implementation prespective we must note,
that in both singlet :eq:`singlet_gamma` and
triplet :eq:`triplet_gamma` a density and magnetic
:math:`\Phi` term appears twice.
Once without an index flip and a dependence on :math:`K-K'`,
:math:`\Phi_{a\overline{b}c\overline{d}}(K-K')`,
and another time with an index flip and a dependence on :math:`K+K'`, 
:math:`\Phi_{c\overline{b}a\overline{d}}(K+K')`.
Inside the linearized Eliashberg equation :eq:`linearized_eliashberg_4`
the :math:`\Phi_{c\overline{b}a\overline{d}}(K+K')` term
picks up a sign which depends on the frequency, momentum and orbital
symmetry of the gap :math:`\Delta`.
For a singlet gap this yields a plus and for a triplet a minus.
Therefore the terms just add up and we lose the factor :math:`1/2` in front of
the square brackets in Eq. :eq:`singlet_gamma` and :eq:`triplet_gamma`.

Also note, that the RPA particle-particle vertices in
Eq. :eq:`singlet_gamma` and :eq:`triplet_gamma` only depend on the difference
between the two fermionic Matsubara frequencies and momenta.
We can therefore write the linearized Eliashberg equation
:eq:`linearized_eliashberg_4` as

.. math::
    \lambda\Delta^{\mathrm{s/t}}(K) =  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}(K-K')
    G(K')G(-K')
    \Delta^{\mathrm{s/t}}(K')\,.
    :label: linearized_eliashberg_5
    
This allows us to get rid of the summation by using the convolution theorem

.. math::
    \lambda
    \mathcal{F}\left[\Delta^{\mathrm{s/t}}(K)\right]=  -\frac{1}{2}
    \mathcal{F}\left[\Gamma^{\mathrm{s/t}}(K-K')\right]
    \mathcal{F}\left[
    G(K')G(-K')
    \Delta^{\mathrm{s/t}}(K')
    \right]\,,
    :label: linearized_eliashberg_5

making the calculation computationaly more efficient.

.. rubric:: References

.. [#abrikosov] A.A. Abrikosov, L.P. Gor’kov, et.al., Pergamon, Oxford (1965)

.. [#yanase] Yanase, et. al., Physics Reports 387, 1-149 (2003)

.. [#takimoto] Takimoto, et. al., PRB 69, 104504 (2004)

.. [#bickers] Bickers, N. E. Self-Consistent Many-Body Theory for Condensed Matter Systems. Theoretical Methods for Strongly Correlated Electrons, 237–296. 6, (2006)

.. [#rohringer] Rohringer, G., New routes towards a theoretical treatment of nonlocal electronic correlations, (2013)

.. [#nourafkan] R. Nourafkan, G. Kotliar, and A. M. Tremblay, Physical Review Letters 117, 1, (Supplementary), (2016) 


