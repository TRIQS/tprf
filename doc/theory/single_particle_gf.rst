.. _single_particle_gf:

On the single particle Green's function
=======================================

The imaginary time single particle Green's function is defined as

.. math::
   G_{a\bar{b}}(\tau_a, \tau_{\bar{b}}) 
   \equiv 
   - \langle \mathcal{T} c_{a}(\tau_a) c^\dagger_{\bar{b}}(\tau_{\bar{b}}) \rangle
   \, .

It is time translational invariant and hence only depends on the time difference

.. math::
   G_{a\bar{b}}(\tau_a, \tau_{\bar{b}}) 
   =
   G_{a\bar{b}}(\tau_a - \tau_{\bar{b}})
   \equiv
   G_{a\bar{b}}(\tau)
   \, .

Using the cyclicity of the trace (see the section on (anti-)periodicity), we can show that for :math:`0 < \tau < \beta`, the bosonic (fermionic) Green's function is :math:`\beta` (anti-)periodic, that is

.. math::
   G_{a\bar{b}}(- \tau) 
   =
   \xi G_{a\bar{b}}(\beta - \tau)

with :math:`\xi = \pm 1` for bosons (fermions).
Hence, extending the function as an (anti-)periodic function to all real valued imaginary times :math:`\tau \in (-\infty, \infty)`, the Green's function can be expanded in the Matsubara Fourier series

.. math::
   G_{a\bar{b}}(\tau) =
   \frac{1}{\beta} \sum_{n=-\infty}^\infty
   e^{- i\nu_n \tau} G_{a\bar{b}}(\nu_n)
   \, ,

with Fourier coefficients

.. math::
   G_{a\bar{b}}(\nu_n) = \int_0^\beta d\tau e^{i\nu_n \tau} G_{a\bar{b}}(\tau)

where :math:`\nu_n` are Matsubara frequencies

.. math::
   \nu_n = \frac{\pi}{\beta}(2n + \vartheta)
   
with :math:`\vartheta = (1-\xi)/2`. From now on, we employ the :math:`\nu \ (\omega)` symbol to denote fermionic (bosonic) Matsubara frequencies.


Field operator Matsubara transforms
-----------------------------------

The notion of the Fourier series can be generalized to the second quantized (field) operators :math:`c(\tau)` and :math:`c^\dagger(\tau)` by introducing the transform relations

.. math::
   c(\nu_n) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{i\nu_n \tau} c(\tau)
   \, , \quad
   c^\dagger(\nu_n) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{-i\nu_n \tau} c^\dagger(\tau)

.. math::
   c(\tau) = \frac{1}{\sqrt{\beta}} \sum_{n=-\infty}^{\infty} e^{-i\nu_n \tau} c(\nu_n)
   \, , \quad
   c^\dagger(\tau) = \frac{1}{\sqrt{\beta}} \sum_{n=-\infty}^{\infty} e^{i\nu_n \tau} c^\dagger(\nu_n)

The symmetic definition of the field operator transforms results in trivial relations for the two-frequency single particle Green's function

.. math::
   G(\nu, \nu') & =
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau e^{i\nu\tau}
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau' e^{-i\nu'\tau'}
   G(\tau, \tau')
   \\ & =
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau e^{i\nu\tau}
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau' e^{-i\nu'\tau'}
   \frac{1}{\beta} \sum_{n=-\infty}^\infty e^{-i \nu''_n (\tau - \tau')}
   G(\nu''_n)
   \\ & =
   \frac{1}{\beta^2} \sum_{n=-\infty}^\infty
   G(\nu''_n)
   \int_0^\beta d\tau e^{(i\nu - i\nu''_n)\tau}
   \int_0^\beta d\tau' e^{(-i\nu' + i\nu''_n)\tau'}
   \\ & =
   \frac{1}{\beta^2} \sum_{n=-\infty}^\infty
   G(\nu''_n)
   \cdot \beta \delta_{\nu, \nu''_n}
   \cdot \beta \delta_{\nu', \nu''_n}
   \\ & =
   \delta_{\nu, \nu'} G(\nu)

Thus there is no scale factor relating the one and two frequency single particle Green's function

.. math::

   G(\nu, \nu') = \delta_{\nu, \nu'} G(\nu)
   \, .
