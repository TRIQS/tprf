.. _boundary_conditions:

(Anti-)Periodicity
==================

Note that the Heisenberg representation the imaginary time dependence of creation and annihilation operators are not conjugated, i.e. :math:`c(\tau) \equiv e^{\tau H} c e^{-\tau H}` and :math:`c^\dagger(\tau) \equiv e^{\tau H} c^\dagger e^{-\tau H}`

.. math::

   G_{a\bar{b}}(\tau)
   \equiv
   - \langle \mathcal{T} a(\tau) \bar{b}(0) \rangle
   \equiv
   - \frac{1}{\mathcal{Z}} \textrm{Tr} \big[
   \mathcal{T} e^{-\int_0^\beta d\bar{\tau} \, H(\bar{\tau})} a(\tau) \bar{b}
   \big]

To derive the boundary conditions we consider two cases. First, for :math:`\beta > \tau > 0` we have

.. math::

   G_{a\bar{b}}(\tau) \Big|_{\beta > \tau > 0}
   =
   - \langle a(\tau) \bar{b}(0) \rangle
   =
   - \frac{1}{\mathcal{Z}}
   \textrm{Tr} \big[
   e^{-\beta H}
   e^{\tau H} a e^{-\tau H}
   \bar{b}
   \big]
   \\ =
   - \frac{1}{\mathcal{Z}}
   \textrm{Tr} \big[
   e^{-\beta H}
   \bar{b}(0)
   e^{(\tau - \beta) H} a e^{-(\tau - \beta) H}
   \big]
   =
   - \langle \bar{b}(0) a(\tau - \beta) \rangle
   \\ =
   - \xi \langle \mathcal{T} a(\tau - \beta) \bar{b}(0) \rangle
     \Big|_{\beta > \tau > 0}
   = \xi G_{a\bar{b}}(\tau - \beta)
     \Big|_{\beta > \tau > 0}

while for :math:`0 > \tau > - \beta` one get

.. math::

   G_{a\bar{b}}(\tau) \Big|_{0 > \tau > -\beta}
   =
   - \xi \langle \bar{b}(0) a(\tau) \rangle
   =
   - \frac{\xi}{\mathcal{Z}}
   \textrm{Tr} \big[
   e^{-\beta H}
   \bar{b}
   e^{\tau H} a e^{-\tau H}
   \big]
   \\ =
   - \frac{\xi}{\mathcal{Z}}
   \textrm{Tr} \big[
   e^{-\beta H}
   e^{(\tau + \beta) H} a e^{-(\tau + \beta) H}
   \bar{b}(0)
   \big]
   =
   - \xi \langle a(\tau + \beta) \bar{b}(0) \rangle
   \\ =
   - \xi \langle \mathcal{T} a(\tau + \beta) \bar{b}(0) \rangle
     \Big|_{0 > \tau > -\beta}
   = \xi G_{a\bar{b}}(\tau + \beta)
     \Big|_{0 > \tau > -\beta}

Thus we see that the single-particle Green's function :math:`G_{a\bar{b}}(\tau)` is :math:`\beta` (anti-)periodic on :math:`\tau \in [\beta, -\beta]`.

.. math::
   G_{a\bar{b}}(\tau) \Big|_{\beta > \tau > 0}
   = \xi G_{a\bar{b}}(\tau - \beta)
     \Big|_{\beta > \tau > 0}   
   \\
   G_{a\bar{b}}(\tau) \Big|_{0 > \tau > -\beta}
   = \xi G_{a\bar{b}}(\tau + \beta)
     \Big|_{0 > \tau > -\beta}

Two-particle Green's functions
------------------------------

The (anti-)periodicity properties can be generalized to two-particle Green's functions in imaginary time. Consider now :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}},\tau_b, \tau_{\bar{c}}, \tau_d=0)`, the cyclic property of the trace and the time-ordering operator (assuming that all operators are either fermionic or bosonic) then yield in the same way

.. math::

   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}}, \tau_d=0)
   \equiv
   \langle \mathcal{T}
   \bar{a}(\tau_{\bar{a}}) b(\tau_b) \bar{c}(\tau_{\bar{c}}) d(0) \rangle

As an example we take the case :math:`\beta > \tau_{\bar{a}} > \tau_b, \tau_{\bar{c}} > 0`

.. math::
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}}, \tau_d=0)
   \equiv
   \langle 
   \bar{a}(\tau_{\bar{a}})
   \big[ \mathcal{T} b(\tau_b) \bar{c}(\tau_{\bar{c}}) \big]
   d(0)
   \rangle
   \\ =
   \langle
   \big[ \mathcal{T} b(\tau_b) \bar{c}(\tau_{\bar{c}}) \big]
   d(0)
   \bar{a}(\tau_{\bar{a}} - \beta)
   \rangle
   =
   \xi \langle
   \mathcal{T}
   \bar{a}(\tau_{\bar{a}} - \beta)
   b(\tau_b) \bar{c}(\tau_{\bar{c}})
   d(0)
   \rangle
   \\ = 
   \xi G^{(2)}_{\bar{a}b\bar{c}d}( \tau_{\bar{a}} - \beta, \tau_b, \tau_{\bar{c}}, 0)
   
In the same way the three periodicity relations read

.. math::
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}})
   \Big|_{\beta > \tau_{\bar{a}} > \tau_b, \tau_{\bar{c}} > 0}
   =
   \xi G^{(2)}_{\bar{a}b\bar{c}d}( \tau_{\bar{a}} - \beta, \tau_b, \tau_{\bar{c}}, 0)
   \\ 
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}})
   \Big|_{\beta > \tau_{b} > \tau_{\bar{a}}, \tau_{\bar{c}} > 0}
   =
   \xi G^{(2)}_{\bar{a}b\bar{c}d}( \tau_{\bar{a}}, \tau_b - \beta, \tau_{\bar{c}}, 0)
   \\ 
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}})
   \Big|_{\beta > \tau_{\bar{c}} > \tau_{\bar{a}}, \tau_{b} > 0}
   =
   \xi G^{(2)}_{\bar{a}b\bar{c}d}( \tau_{\bar{a}}, \tau_b, \tau_{\bar{c}} - \beta, 0)

and the second triple of relations become

.. math::
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}})
   \Big|_{\beta > \tau_b, \tau_{\bar{c}} > 0 > \tau_{\bar{a}} > -\beta}
   =
   \xi G^{(2)}_{\bar{a}b\bar{c}d}( \tau_{\bar{a}} + \beta, \tau_b, \tau_{\bar{c}}, 0)
   \\ 
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}})
   \Big|_{\beta > \tau_{\bar{a}}, \tau_{\bar{c}} > 0 > \tau_{b} > -\beta}
   =
   \xi G^{(2)}_{\bar{a}b\bar{c}d}( \tau_{\bar{a}}, \tau_b + \beta, \tau_{\bar{c}}, 0)
   \\ 
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}})
   \Big|_{\beta > \tau_{\bar{a}}, \tau_{b} > 0 > \tau_{\bar{c}} > -\beta}
   =
   \xi G^{(2)}_{\bar{a}b\bar{c}d}( \tau_{\bar{a}}, \tau_b, \tau_{\bar{c}} + \beta, 0)
   
Kubo-Martin-Schwinger (KMS) boundary conditions
===========================================================

The boundary conditions in imaginary time for the Green's functions are generated by  the commutation relation :math:`[a, \bar{b}]_{-\xi} = a\bar{b} - \xi \bar{b}a = \delta_{ab}`, where :math:`\xi = \pm 1` for bosons and fermions respectively

.. math::
   G_{a\bar{b}}(0^+) = -\langle \mathcal{T} a(0^+) \bar{b}(0) \rangle
   = -\langle a \bar{b} \rangle
   \\
   G_{a\bar{b}}(0^-) = -\langle \mathcal{T} a(0^-) \bar{b}(0) \rangle
   = - \xi \langle \bar{b} a \rangle
   
So that the boundary condition at :math:`\tau = 0^\pm` is

.. math::
   G_{a\bar{b}}(0^+) - G_{a\bar{b}}(0^-)
   = -\langle a \bar{b} - \xi \bar{b} a \rangle
   = -\langle [a, \bar{b}]_{-\xi} \rangle
   = -\delta_{ab}

Using the periodicity relation :math:`G_{a\bar{b}}(0^-) = \xi G_{a\bar{b}}(\beta^-)` we finally arrive at the boundary condition restricted to :math:`\beta > \tau > 0`
   
.. math::
   G_{a\bar{b}}(0^+) - \xi G_{a\bar{b}}(\beta^-) = - \delta_{ab}

.. note::
   The anomalous Green's functions has the simpler boundary condition

   .. math::
      G_{ab}(0^+) - \xi G_{ab}(\beta^-) = 0
      \\
      G_{\bar{a}\bar{b}}(0^+) - \xi G_{\bar{a}\bar{b}}(\beta^-) = 0

   since :math:`[a, b]_{-\xi} = 0` and :math:`[\bar{a}, \bar{b}]_{-\xi} = 0`.

Generalized Kubo-Martin-Schwinger (KMS) boundary conditions
===========================================================

For the two-particle Green's function the KMS boundary conditions generalize to relations incorporating the single particle Green's function.

.. math::
   G^{(2)}_{\bar{a}b\bar{c}d}(0^+,\tau_b, \tau_{\bar{c}})
   - \xi
   G^{(2)}_{\bar{a}b\bar{c}d}(\beta^-,\tau_b, \tau_{\bar{c}})
   = \xi \delta_{ad} G_{b\bar{c}}(\tau_b - \tau_{\bar{c}})
   \\
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, 0^+)
   - \xi
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau_b, \beta^-)
   = \delta_{cd} G_{b\bar{a}}(\tau_b - \tau_{\bar{a}})
   \\
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, 0^+, \tau_{\bar{c}})
   - \xi
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \beta^-, \tau_{\bar{c}})
   = 0

Thus the discontinuities at :math:`\tau_{\bar{a}}=0` and :math:`\tau_{\bar{c}}=0` are non-trivial and given by the single-particle Green's function.

The two additional discontinuities in :math:`\tau_{\bar{a}}, \tau_b, \tau_{\bar{c}} \in [\beta, 0]` are the three equal time planes :math:`\tau_{\bar{a}} = \tau_b`, :math:`\tau_b = \tau_{\bar{c}}`, and :math:`\tau_{\bar{a}} = \tau_{\bar{c}}`.

.. math::

   G^{(2)}_{\bar{a}b\bar{c}d}(\tau^+, \tau^-, \tau_{\bar{c}})
   -
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau^-, \tau^+, \tau_{\bar{c}})
   = -\delta_{bc} G_{d\bar{a}}(\beta - \tau_{\bar{a}})
   \\
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau^+, \tau^-)
   -
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau_{\bar{a}}, \tau^-, \tau^+)
   = - \delta_{bc} G_{d\bar{a}}( \beta - \tau_{\bar{a}} )
   \\
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau^+, \tau_b, \tau^-)
   -
   G^{(2)}_{\bar{a}b\bar{c}d}(\tau^-, \tau_b, \tau^+)
   = 0

.. note::
   Write test on `pyed` output to explicitly check the boundary conditions.

   Figure out how to enforce boundary conditions on stochastically sampled data.
   
   
