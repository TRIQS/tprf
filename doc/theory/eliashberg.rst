.. _eliashberg:

Linearized Eliashberg Equation
================================

The linearized Eliashberg equation is a generalization of the linearized
Bardeen-Cooper-Schrieffer (BCS) gap equation to frequency dependent gaps.
It can be used to determine the critical temperature :math:`T_\mathrm{c}`,
at which a transition to a superconducting state occurs,
and the symmetry of the corresponding gap function :math:`\Delta`.
It is given by

.. math::
    \Delta_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{PP}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
    \Delta_{\bar{e}\bar{f}}(K')\,.
    :label: linearized_eliashberg_1

where :math:`Q/K` is a combination of bosonic/fermionic Matsubara 
frequency and momentum,
:math:`N_{\mathbf{k}}` is the number of momentum points,
:math:`\Gamma^{PP}` is the irreducible particle-particle vertex
and :math:`G` is the one-particle Green's function.

Note, that the bosonic Matsubara frequency and momentum in the particle-particle vertex
is set to zero.
This is because we are interested in Cooper-pairs which have a zero
transfered momentum-frequency in a scattering process.
    
Deriving the linearized Eliashberg equation: Normal state and superconducting state
-----------------------------------------------------------------------------------

The linearized Eliashberg equation can be seen from two perspectives.
On one hand we are in the normal state and want to find the transition
to the superconducting one,
and on the other we are in the superconducting state using the limit of
small gaps :math:`\Delta \ll 1`.

Normal state
^^^^^^^^^^^^

Generally speaking a transition from the normal state to a superconducting
one occurs when the particle-particle susceptibility diverges.

.. math::
    \mathbf{\chi}^{PP} = [\mathbf{1}-\mathbf{\Gamma}^{PP} \mathbf{\chi}^{(0),{PP}}]^{-1}
    \mathbf{\chi}^{(0),{PP}}

This is the case when the largest eigenvalue of 
:math:`\mathbf{\Gamma}^{PP} \mathbf{\chi}^{(0),{PP}}` becomes unity.
For a largest eigenvalues that is smaller than :math:`1` we are still in the
normal state,
but we can calculate the corresponding eigenvectors :math:`\Delta`.
This corresponds to the following eigenvalue equation

.. math::
    \lambda\Delta_{\bar{a}\bar{b}}(K)= \frac{T^2_{\mathrm{c}}}{2 N_{\mathbf{k}}^2}\sum_{K', K''}
    \Gamma^{PP}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    \chi^{(0),{PP}}_{\bar{e}d\bar{f}c}(Q=0, K', K'')
    \Delta_{\bar{e}\bar{f}}(K')\,,
    :label: linearized_eliashberg_2

which we can write like Eq. :eq:`linearized_eliashberg_1` with the definiton
of :math:`\chi^{(0),{PP}}` :eq:`bare_pp_sus_def`, 

.. math::
    \lambda\Delta_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{PP}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
    \Delta_{\bar{e}\bar{f}}(K')\,.

This equation is valid for :math:`\lambda \leq 1`
and yields eigenvectors, which are superconducting gaps that have not manifested yet.
At :math:`\lambda=1` the normal state breaks down and the superconducting
state with the corresponding gap emerges.
The size of eigenvalues is therefore an indicator of how likely the associated gap is
to manifest.

Superconducting state
^^^^^^^^^^^^^^^^^^^^^

For a calculation in the superconducting state we must extend the basis
to accommodate Cooper-pairs.
In addition to the normal one-particle Green's function

.. math::
   G_{a\bar{b}}(\tau, \mathbf{r}) 
   \equiv 
   - \langle \mathcal{T} c_{a\uparrow}(\tau, \mathbf{r})
   c^\dagger_{\bar{b}\uparrow}(0, \mathbf{0}) \rangle \,,

and its backwards propagating counterpart

.. math::
   \overline{G}_{\bar{a}b}(\tau, \mathbf{r}) 
   \equiv 
   - \langle \mathcal{T} c^\dagger_{\bar{a}\downarrow}(\tau, \mathbf{r})
   c_{b\downarrow}(0, \mathbf{0}) \rangle\,,

we have to introduce the one-particle anomalous Green's functions
:math:`F` and :math:`\overline{F}`.
These are defined as

.. math::
    F_{ab}(\tau, \mathbf{r}) 
    \equiv
   \langle \mathcal{T} c_{a\uparrow}(\tau, \mathbf{r}) 
   c_{b\downarrow}(0, \mathbf{0}) \rangle
   \,,

and

.. math::
   \overline{F}_{\bar{a}\bar{b}}(\tau, \mathbf{r}) 
   \equiv
   \langle \mathcal{T} c^\dagger_{\bar{a}\downarrow}(\tau, \mathbf{r}) 
   c^\dagger_{\bar{b}\uparrow}(0, \mathbf{0}) \rangle\,.

Fourier transforming to Matsubara frequency and momentum space then gives that

.. math::
   \overline{G}_{\bar{a}b}(i\nu_n, \mathbf{k}) = 
   -G_{b\bar{a}}(-i\nu_n, -\mathbf{k})\,,
   :label: g_bar_to_g

and

.. math::
   \overline{F}_{\bar{a}\bar{b}}(i\nu_n, \mathbf{k})
   =
   [F_{ab}(i\nu_n, \mathbf{k}) ]^{\dagger}\,.
   :label: f_bar_to_f

All four Green's functions are coupled and given by

.. math::
   \left( \begin{array}{cc}
     \mathbf{G} & \mathbf{F} \\
     \overline{\mathbf{F}} & \overline{\mathbf{G}}\\
   \end{array} \right)
   =
   \left( \begin{array}{cc}
     \left(\mathbf{G}^{(0)}\right)^{-1} - \mathbf{\Sigma} & \mathbf{\Delta} \\
     \overline{\mathbf{\Delta}} & \left(\overline{\mathbf{G}}^{(0)}\right)^{-1} - \overline{\mathbf{\Sigma}} \\
   \end{array} \right)^{-1}\,,

where :math:`\Sigma/\overline{\Sigma}` are the normal self-energies
and :math:`\Delta/\overline{\Delta}` the anomalous ones.
The anomalous self-energies are equivalent to the superconducting gap and are given by

.. math::
    \Delta_{\bar{a}\bar{b}}(K)
    =
    \frac{1}{2N_{\mathbf{k}} \beta} \sum_{K'}
    \Gamma^{PP}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    F_{cd}(K')\,,
    :label: anomalous_self_energy

.. math::
    \overline{\Delta}_{{a}{b}}(K)
    =
    \frac{1}{2N_{\mathbf{k}} \beta} \sum_{K'}
    \Gamma^{PP}_{a\bar{c}b\bar{d}}(Q=0, K, K')
    \overline{F}_{\bar{c}\bar{d}}(K')\,.
    :label: anomalous_self_energy_2

With either of those equations we could calculate the gap for the superconducting state
below :math:`T_\mathrm{c}`.
But note, that this is not trivial, due to the coupling of the Green's functions and
self-energies.
In the limit close to :math:`T_\mathrm{c}` the superconducting gap is very small and
we can approximate the anomalous Green's function with

.. math::
    \mathbf{F} = -
    \left[\left(\mathbf{G}^{(0)}\right)^{-1} - \mathbf{\Sigma}\right]^{-1}
    \mathbf{\Delta}
    \left[\left(\mathbf{\overline{G}}^{(0)}\right)^{-1} - \mathbf{\overline{\Sigma}}\right]^{-1}
    :label: approx_anomalous_gf

Plugging this linearized version of the anomalouse Green's function in
Eq. :eq:`anomalous_self_energy` with the help of relation
Eq. :eq:`g_bar_to_g` yields the linearized Eliashberg equation
Eq. :eq:`linearized_eliashberg_1`.

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
    \Delta(\mathbf{k}) =  -\frac{1}{2 N_{\mathbf{k}}}\sum_{\mathbf{k'}}
    \Gamma^{PP}(\mathbf{q}=\mathbf{0}, \mathbf{k}, \mathbf{k'})
    \frac{\tan(\epsilon(\mathbf{k'})\beta/2)}{2\epsilon(\mathbf{k'})}
    \Delta(\mathbf{k'})\,,
    :label: linearized_eliashberg_3

which corresponds to the linearized BCS gap equation.
The non-linear BCS gap equation can be obtained from Eq. :eq:`linearized_eliashberg_3` 
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
it can be either even (:math:`+`) or odd (:math:`-`) under each separate operation,
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
under these four degrees of freedom. We list all eight possible combinations
in the table below: 

.. table:: 
    :align: center
    :widths: grid

    +-----------------+-----------------+-----------------+-----------------+
    |        S        |        P        |        O        |        T        |
    +=================+=================+=================+=================+
    |    :math:`-`    |    :math:`+`    |    :math:`+`    |    :math:`+`    |
    +-----------------+-----------------+-----------------+-----------------+
    |    :math:`-`    |    :math:`-`    |    :math:`-`    |    :math:`+`    |
    +-----------------+-----------------+-----------------+-----------------+
    |    :math:`-`    |    :math:`-`    |    :math:`+`    |    :math:`-`    |
    +-----------------+-----------------+-----------------+-----------------+
    |    :math:`-`    |    :math:`+`    |    :math:`-`    |    :math:`-`    |
    +-----------------+-----------------+-----------------+-----------------+
    |    :math:`+`    |    :math:`-`    |    :math:`-`    |    :math:`-`    |
    +-----------------+-----------------+-----------------+-----------------+
    |    :math:`+`    |    :math:`+`    |    :math:`+`    |    :math:`-`    |
    +-----------------+-----------------+-----------------+-----------------+
    |    :math:`+`    |    :math:`+`    |    :math:`-`    |    :math:`+`    |
    +-----------------+-----------------+-----------------+-----------------+
    |    :math:`+`    |    :math:`-`    |    :math:`+`    |    :math:`+`    |
    +-----------------+-----------------+-----------------+-----------------+

Because all other combinations are unphysical it is possible to restrict the gap to the
allowed symmetries while solving the linearized Eliashberg equation. 


Spin diagonalization
^^^^^^^^^^^^^^^^^^^^

For :math:`SU(2)` symmetric systems we can drop the spin dependency
by diagonalizing everything in spin.
This diagonalization splits the superconducting gap :math:`\Delta`
in two channels,
the singlet channel

.. math::
    \Delta^{\mathrm{s}}
    =
    \Delta_{\uparrow\downarrow}
    -
    \Delta_{\downarrow\uparrow}\,,

and triplet channel

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

This is analog for the triplet channel and we obtain the spin diagonalized
linearized Eliashberg equation

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(Q=0, K, K')
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,,
    :label: linearized_eliashberg_4

with all indices being only orbital ones.

Random phase approximation for the spin diagonalized irreducible particle-particle vertex
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To obtain the spin diagonalized irreducible particle-particle vertex in the
random phase approximation (RPA)
one substitutes all vertices with the bare one in the parquet equation.
This yields for the singlet channel

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

and for the triplet channel

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
	\Phi^{\mathrm{d/m}}_{a\bar{b}c\bar{d}}
	\equiv
	U_{a\bar{b}e\bar{f}}^{\mathrm{d/m}}
    \chi_{\bar{f}e\bar{g}h}^{\mathrm{d/m}}(Q)
    U^{\mathrm{d/m}}_{h\bar{g}c\bar{d}}
    \,.

Note, that the superscripts :math:`\mathrm{d}` and :math:`\mathrm{m}`
indicate the density and magnetic channel.
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
Inside the linearized Eliashberg equation :eq:`linearized_eliashberg_4`
the :math:`\Phi_{c\overline{b}a\overline{d}}(K+K')` term
picks up a sign which depends on the frequency, momentum and orbital
symmetry of the gap :math:`\Delta`.
For all allowed singlet combinations it is positive and for all allowed triplet ones
negative. Therefore Eq. :eq:`singlet_gamma` and Eq. :eq:`triplet_gamma` become

.. math::
    \Gamma^{\text{s}}_{a\overline{b}c\overline{d}}(Q=0, K, K') \equiv
    \Lambda^{\text{s}}_{a\overline{b}c\overline{d}}
    +
    3 
    \Phi^{\text{m}}_{a\overline{b}c\overline{d}}(K-K')
    -
    \Phi^{\text{d}}_{a\overline{b}c\overline{d}}(K-K')
    \,,
    :label: singlet_gamma_2

.. math::
    \Gamma^{\text{t}}_{a\overline{b}c\overline{d}}(Q=0, K, K') \equiv
    \Lambda^{\text{t}}_{a\overline{b}c\overline{d}} 
    -
    \Phi^{\text{m}}_{a\overline{b}c\overline{d}}(K-K')
    -
    \Phi^{\text{d}}_{a\overline{b}c\overline{d}}(K-K')
    \,.
    :label: triplet_gamma_2

Note, that this simplification is only allowed, if the solutions of :math:`\Delta`
are restricted to the allowed symmetries, otherwise unphysical solution can occur.
Also note, that the RPA particle-particle vertices in
Eq. :eq:`singlet_gamma_2` and :eq:`triplet_gamma_2` only depend on the difference
between the two fermionic Matsubara frequencies and momenta.
We can therefore write the linearized Eliashberg equation
:eq:`linearized_eliashberg_4` as

.. math::
    \lambda\Delta^{\mathrm{s/t}}_{\bar{a}\bar{b}}(K)=  -\frac{1}{2 N_{\mathbf{k}}\beta}\sum_{K'}
    \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(K-K')
    G_{c\bar{e}}(K')G_{d\bar{f}}(-K')
    \Delta^{\mathrm{s/t}}_{\bar{e}\bar{f}}(K')\,.
    :label: linearized_eliashberg_5
    
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
