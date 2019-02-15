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

With the arise of Cooper pairs we need in addition to the normal single-particle Green's function

.. math::
   G_{a\bar{b}}(\tau_1, \tau_2) 
   \equiv 
   - \langle \mathcal{T} c_{a}(\tau_1) c^\dagger_{\bar{b}}(\tau_2) \rangle\,,

and its backwards propagating counterpart

.. math::
   \bar{G}_{\bar{a}b}(\tau_1, \tau_2) 
   \equiv 
   - \langle \mathcal{T} c^\dagger_{\bar{a}}(\tau_1) c_{b}(\tau_2) \rangle\,,

the single-particle anomalous Green's functions :math:`F` and :math:`\bar{F}` to describe a superconducting state.
These are defined as

.. math::
    F_{a\bar{b}}(\tau_1, \tau_2) 
    \equiv
   \langle \mathcal{T} c_{a}(\tau_1) c_{\bar{b}}(\tau_2) \rangle\,,

.. math::
    \bar{F}_{a\bar{b}}(\tau_1, \tau_2) 
    \equiv
   \langle \mathcal{T} c^\dagger_{a}(\tau_1) c^\dagger_{\bar{b}}(\tau_2) \rangle\,.

Dyson-Gorkov Equations
----------------------

The former properties of a superconductor are given by four coupled equations called Dyson-Gorkov Equations.

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
The anomalous self-energy can be expressed with the effective pairing interaction :math:`\Gamma` and the anomalous Green's function :math:`F` as

.. math::
    \Delta_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{q}} \sum_{i\nu'}
    \Gamma_{A\bar{a}\bar{b}B}(\mathbf{k}-\mathbf{q}, i\nu - i\nu') F_{AB}(\mathbf{q}, i\nu')\,.
    :label: anom_self_energy

Around the transition point to the superconducting state the anomalous self-energy :math:`\Delta` is approximately zero, and, because we are only interested in the transition point, we linearize the Dyson-Gorkov equations with respect to :math:`\Delta`.
This yields

.. math::
    G(a\bar{b}) = \big({G^{(0)}}^{-1}(a\bar{b}) - \Sigma(a\bar{b})\big)^{-1}\,,

.. math::
    \bar{G}(\bar{a}b) = \big({\bar{G}^{(0)}}^{-1}(\bar{a}b) - \bar{\Sigma}(\bar{a}b)\big)^{-1}\,,

.. math:: 
    F(ab) =  \big({G^{(0)}}^{-1}(a\bar{c}) - \Sigma(a\bar{c})\big)^{-1} \Delta(\bar{c}\bar{d}) 
        \big({\bar{G}^{(0)}}^{-1}(\bar{d}b) - \bar{\Sigma}(\bar{d}b)\big)^{-1}\,,
    :label: lin_anom_gf

.. math:: 
    \bar{F}(\bar{a}\bar{b}) =  \big({\bar{G}^{(0)}}^{-1}(\bar{a}c) - \bar{\Sigma}(\bar{a}c)\big)^{-1} \Delta(cd) 
        \big({G^{(0)}}^{-1}(d\bar{b}) - \Sigma(d\bar{b})\big)^{-1}\,.

We can then plug :eq:`lin_anom_gf` into :eq:`anom_self_energy` to obtain the linearized Eliashberg equation

.. math::
    \Delta_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{q}} \sum_{i\nu'}
    \Gamma_{A\bar{a}\bar{b}B}(\mathbf{k}-\mathbf{q}, i\nu - i\nu') 
    \big({G^{(0)}}^{-1}(\mathbf{q}, i\nu') - \Sigma(\mathbf{q}, i\nu') \big)^{-1}_{A\bar{c}} \\
    \Delta_{\bar{c}\bar{d}}(\mathbf{q}, i\nu')
    \big({G^{(0)}}^{-1}_{}(-\mathbf{q}, -i\nu') - \Sigma_{}(-\mathbf{q}, -i\nu') \big)^{-1}_{B\bar{d}}\,.

To make use of this equations it is usually interpreted as an eigenvalue equation

.. math::
    \lambda \Delta = \Lambda \Delta\,,

where the eigenvalue :math:`\lambda` is seen as a measurement for the strength of superconducting ordering and a phase transition occurs when it reaches unity.

RPA Approach
------------

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





