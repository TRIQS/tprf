.. _eliashberg:

Linearized Eliashberg Equation
================================

.. note::
    References:

    - [A.A. Abrikosov, L.P. Gor’kov, et.al., Pergamon, Oxford (1965)]
    - [Takimoto, et. al., PRB 69, 104504 (2004)]
    - [Yanase, et. al., Physics Reports 387, 1-149 (2003)]


.. note::
    All indices on this page only represent orbital degrees of freedom.
    Spin is not treated explicitly and therefore only spin-independent Hamiltonians can be used for calculations. 




We asssume a homogenous system with some arbitrary effective pairing interaction :math:`\Gamma`, which leads to the formation of Cooper pairs.

Anomalous Green's Functions
---------------------------

.. note::
   Explain what happens with all spin quantum numbers in the single-particle Green's function. Do we work with a particular combination of spins? :math:`G_{a\bar{b}} = G_{\alpha \uparrow \bar{b} \downarrow}`?

With the arise of Cooper pairs we need in addition to the normal single-particle Green's function

.. math::
   G_{a\bar{b}}(\tau - \tau') 
   \equiv 
   - \langle \mathcal{T} c_{a}(\tau) c^\dagger_{\bar{b}}(\tau') \rangle
   =
   - \langle \mathcal{T} a(\tau) \bar{b}(\tau') \rangle\,,

and its backwards propagating counterpart

.. math::
   \bar{G}_{\bar{a}b}(\tau - \tau') 
   \equiv 
   - \langle \mathcal{T} c^\dagger_{\bar{a}}(\tau) c_{b}(\tau') \rangle
   =
   - \langle \mathcal{T} \bar{a}(\tau) b(\tau') \rangle\,,

the single-particle anomalous Green's functions :math:`F` and :math:`\bar{F}` to describe a superconducting state.
These are defined as

.. math::
    F_{ab}(\tau - \tau') 
    \equiv
   \langle \mathcal{T} c_{a}(\tau) c_{b}(\tau') \rangle
   =
   \langle \mathcal{T} a(\tau) b(\tau') \rangle
   \,,

.. math::
   \bar{F}_{\bar{a}\bar{b}}(\tau - \tau') 
   \equiv
   \langle \mathcal{T} c^\dagger_{a}(\tau) c^\dagger_{\bar{b}}(\tau') \rangle
   =
   \langle \mathcal{T} \bar{a}(\tau) \bar{b}(\tau') \rangle\,.

Fourier transforming to Matsubara frequency space then gives that

.. math::
   \bar{G}_{\bar{a}b}(\mathbf{k}, i\nu_n) = [ G_{b\bar{a}}(-\mathbf{k}, -i\nu_n) ]^{*}
   \\
   \bar{F}_{\bar{a}\bar{b}}(\mathbf{k}, i\nu_n) = [ F_{ba}(-\mathbf{k}, -i\nu_n) ]^{*} 

Dyson-Gorkov Equations
----------------------

The former properties of a superconductor are given by the Dyson-Gorkov equations

.. math::
   \mathbf{G}(\mathbf{k}, i\nu_n)
   =
   \mathbf{G}^{(0)}(\mathbf{k}, i\nu_n)
   + \mathbf{G}^{(0)}(\mathbf{k}, i\nu_n)
     \ast \mathbf{\Sigma}(\mathbf{k}, i\nu_n)
     \ast \mathbf{G}(\mathbf{k}, i\nu_n)

.. math::
   \mathbf{G} \equiv
   \left[ \begin{array}{cc}
     G_{a\bar{b}} & F_{ab} \\
     \bar{F}_{\bar{a}\bar{b}} & \bar{G}_{\bar{a}b} \\
   \end{array} \right]
   \quad
   \mathbf{G}^{(0)}
   \equiv
   \left[ \begin{array}{cc}
     G^{(0)}_{a\bar{b}} & 0 \\
     0 & \bar{G}^{(0)}_{\bar{a}b} \\
   \end{array} \right]
   \quad
   \mathbf{\Sigma}
   \equiv
   \left[ \begin{array}{cc}
     \Sigma_{a\bar{b}} & \Delta_{ab} \\
     \bar{\Delta}_{\bar{a}\bar{b}} & \bar{\Sigma}_{\bar{a}b} \\
   \end{array} \right]

In component form this becomes,

.. math::
    G(a\bar{b}) = G^{(0)}(a\bar{b}) + G^{(0)}(a\bar{c})\Sigma(\bar{c}d)G(d\bar{b}) +
    G^{(0)}(a\bar{c})\bar{\Delta}(\bar{c}\bar{d})\bar{F}(\bar{d}\bar{b})

.. math::
    \bar{G}(\bar{a}b) = \bar{G}^{(0)}(\bar{a}b) + \bar{G}^{(0)}(\bar{a}c)\bar{\Sigma}(c\bar{d})\bar{G}(\bar{d}b) +
    \bar{G}^{(0)}(\bar{a}c)\Delta(cd)F(db)

.. math::
    F(ab) = G^{(0)}(a\bar{c}) \Sigma(\bar{c}d) F(db)+
    G^{(0)}(a\bar{c}) \bar{\Delta}(\bar{c}\bar{d}) \bar{G}(\bar{d}b)

.. math::
    \bar{F}(\bar{a}\bar{b}) = \bar{G}^{(0)}(\bar{a}c) \bar{\Sigma}(c\bar{d}) \bar{F}(\bar{d}\bar{b})+
    \bar{G}^{(0)}(\bar{a}c) \Delta(cd) G(d\bar{b})

Here :math:`\Sigma` is the normal self-energy and :math:`\Delta` and :math:`\bar{\Delta}` the anomalous self-energies, which are equal in the absence of a magnetic field and will be treated as from now on.

Anomalous self-energy and particle-particle vertex
--------------------------------------------------

.. note::
   Define :math:`\Gamma`. It should be the particle-particle vertex :math:`\Gamma^{(pp)}` related to the generalized susceptibility :math:`\chi` through the Bethe-Salpeter equation in the particle-particle channel. This would give the definition of the four orbital(spin) indices and their order.

The anomalous self-energy can be expressed with the effective pairing interaction :math:`\Gamma` and the anomalous Green's function :math:`F` as

.. math::
    \Delta_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{q}} \sum_{i\nu'}
    \Gamma_{A\bar{a}\bar{b}B}(\mathbf{k}-\mathbf{q}, i\nu - i\nu') F_{AB}(\mathbf{q}, i\nu')\,.
    :label: anom_self_energy

Linearization in :math:`\Delta`
-------------------------------
	    
Around the transition point to the superconducting state the anomalous self-energy :math:`\Delta` is approximately zero, and, because we are only interested in the transition point, we linearize :math:`F` in the Dyson-Gorkov equations with respect to :math:`\Delta`. This yields

.. math::
   F & = g \Sigma F + g \Delta \bar{G} \\
   G & = g + g \Sigma G + g \Delta \bar{F}

.. math::
   F & = (g^{-1} - \Sigma)^{-1} \Delta \bar{G} \\
   \bar{G} & = (\bar{g}^{-1} - \Sigma)^{-1} + \bar{\Delta} F

.. math::
   F = (g^{-1} - \Sigma)^{-1} \Delta (\bar{g}^{-1} - \bar{\Sigma})^{-1} + \mathcal{O}(\Delta^2)
   :label: lin_anom_gf
   
We then insert :eq:`lin_anom_gf` into :eq:`anom_self_energy` and obtain the linearized Eliashberg equation

.. math::
    \Delta_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{q}} \sum_{i\nu'}
    \Gamma_{A\bar{a}\bar{b}B}(\mathbf{k}-\mathbf{q}, i\nu - i\nu')
    \\ \times
    \big({G^{(0)}}^{-1}(\mathbf{q}, i\nu') - \Sigma(\mathbf{q}, i\nu') \big)^{-1}_{A\bar{c}}
    \Delta_{\bar{c}\bar{d}}(\mathbf{q}, i\nu')
    \big({G^{(0)}}^{-1}_{}(-\mathbf{q}, -i\nu') - \Sigma_{}(-\mathbf{q}, -i\nu') \big)^{-1}_{B\bar{d}}\,.

To make use of this equations it is usually interpreted as an eigenvalue equation

.. math::
    \lambda \Delta = \Lambda \Delta\,,

where the eigenvalue :math:`\lambda` is seen as a measurement for the strength of superconducting ordering and a phase transition occurs when it reaches unity.

RPA Approach
------------

.. note::
   Explain what happens with momenta

The linearized Eliashberg equation can be studied in the RPA limit.
In this case the normal self-energy is set to zero and the effective pairing interaction :math:`\Gamma` for triplet Cooper pairs is given by

.. math::
    \Gamma^{(\mathrm{triplet})}(\bar{a}b\bar{c}d) =
    \frac{3}{2} U^{(\mathrm{s})}(\bar{a}b\bar{A}B) \chi^{(\mathrm{s})}(\bar{A}B\bar{C}D) 
    U^{(\mathrm{s})}(\bar{C}D\bar{c}d) \\
    -\frac{1}{2} U^{(\mathrm{c})}(\bar{a}b\bar{A}B) \chi^{(\mathrm{c})}(\bar{A}B\bar{C}D) 
    U^{(\mathrm{c})}(\bar{C}D\bar{c}d) \\
   + \frac{1}{2}\big(U^{(\mathrm{s})}(\bar{a}b\bar{c}d)+
    U^{(\mathrm{c})}(\bar{a}b\bar{c}d)\big)\,.

Here :math:`\chi^{(\mathrm{s})}` is the spin-susceptibility tensor defined by

.. math::
    \chi^{(\mathrm{s})}(\bar{a}b\bar{c}d) = \big(\mathbb{1} - \chi^{(0)}(\bar{a}b\bar{A}B) 
    U^{(\mathrm{s})}(\bar{A}B\bar{C}D)\big)^{-1}  \chi^{(0)}(\bar{C}D\bar{c}d)\,,

and :math:`\chi^{(\mathrm{c})}` is the charge-susceptibility tensor defined by

.. math::
    \chi^{(\mathrm{c})}(\bar{a}b\bar{c}d) = \big(\mathbb{1} + \chi^{(0)}(\bar{a}b\bar{A}B) 
    U^{(\mathrm{c})}(\bar{A}B\bar{C}D)\big)^{-1}  \chi^{(0)}(\bar{C}D\bar{c}d)\,,

here :math:`\chi^{(0)}` is the non-interacting particle-hole bubble.

The spin and charge interaction tensors are given by

.. math::
    U^{(\mathrm{s})}(a\bar{a}b\bar{b}) =
    \begin{cases}
    U, & \mathrm{if}\;a=\bar{a}=b=\bar{b} \\
    U', & \mathrm{if}\;a=\bar{b}\neq \bar{a}=b \\
    J, & \mathrm{if}\;a=\bar{a}\neq b=\bar{b} \\
    J', & \mathrm{if}\;a=b\neq \bar{a}=\bar{b} \\
    0, & \mathrm{else}
    \end{cases}


.. math::
    U^{(\mathrm{c})}(a\bar{a}b\bar{b}) =
    \begin{cases}
    U, & \mathrm{if}\;a=\bar{a}=b=\bar{b} \\
    -U'+2J, & \mathrm{if}\;a=\bar{b}\neq \bar{a}=b \\
    2U'-J, & \mathrm{if}\;a=\bar{a}\neq b=\bar{b} \\
    J', & \mathrm{if}\;a=b\neq \bar{a}=\bar{b} \\
    0, & \mathrm{else}
    \end{cases}

where :math:`U`, :math:`U'`, :math:`J` and :math:`J'` are the usual Kanamori interaction parameters.





