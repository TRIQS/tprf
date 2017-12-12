.. _single_particle_gf:

On the single particle Green's function
=======================================

The imaginary time single particle Green's function is defined as

.. math::
   G_{a\bar{b}}(\tau_1, \tau_2) 
   \equiv 
   - \langle \mathcal{T} c_{a}(\tau_1) c^\dagger_{\bar{b}}(\tau_2) \rangle
   \, .

It is time translational invariant and hence only depends on the time difference

.. math::
   G_{a\bar{b}}(\tau_1, \tau_2) 
   =
   G_{a\bar{b}}(\tau_1 - \tau_2)
   \, .

Using the cyclicity of the trace, see the section on (anti-)periodicity, we can show that it is :math:`\beta` (anti-)periodic.

.. math::
   G_{a\bar{b}}(- \tau) 
   =
   \xi G_{a\bar{b}}(\beta - \tau)
   \, .

Hence, extending the function as an (anti-)periodic function to all real valued imaginary times :math:`\tau \in (-\infty, \infty)` the Green's function can be expanded in the Matsubara Fourier series

.. math::
   G_{a\bar{b}}(\tau) =
   \frac{1}{\beta} \sum_{n=-\infty}^\infty
   e^{- i\nu_n \tau} G_{a\bar{b}}(i\nu_n)
   \, ,

with Fourier coefficients

.. math::
   G(i\nu_n) = \int_0^\beta d\tau e^{i\nu_n \tau} G_{a\bar{b}}(\tau)

where :math:`\nu_n` are Matsubara frequencies

.. math::
   \nu_n = \frac{\pi}{\beta}(2n + \vartheta)
   
with :math:`\vartheta = (1-\xi)/2` and :math:`\xi = \pm 1` for bosons/fermions) exploiting :math:`\beta` (anti)periodicity.


Field operator Matsubara transforms
-----------------------------------

The notion of the Fourier series can be generalized to the second quantized (field) operators :math:`c(\tau)` and :math:`c^\dagger(\tau)` by introducing the transform relations

.. math::
   c(i\nu_n) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{i\nu_n \tau} c(\tau)
   \, , \quad
   c^\dagger(i\nu_n) \equiv \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau \, e^{-i\nu_n \tau} c^\dagger(\tau)

.. math::
   c(\tau) = \frac{1}{\sqrt{\beta}} \sum_{n=-\infty}^{\infty} e^{-i\nu_n \tau} c(i\nu_n)
   \, , \quad
   c^\dagger(\tau) = \frac{1}{\sqrt{\beta}} \sum_{n=-\infty}^{\infty} e^{i\nu_n \tau} c^\dagger(i\nu_n)

The symmetic definition of the field operator transforms results in trivial relations for the two frequency single particle Green's function

.. math::
   G(i\nu, i\nu') =
   \\ =
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau e^{i\nu\tau}
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau' e^{-i\nu'\tau'}
   G(\tau, \tau')
   \\ =
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau e^{i\nu\tau}
   \frac{1}{\sqrt{\beta}} \int_0^\beta d\tau' e^{-i\nu'\tau'}
   \frac{1}{\beta} \sum_{n=-\infty}^\infty e^{-i \omega (\tau - \tau')}
   G(i\omega)
   \\ =
   \frac{1}{\beta^2} \sum_{n=-\infty}^\infty
   G(i\omega)
   \int_0^\beta d\tau e^{(i\nu - i\omega)\tau}
   \int_0^\beta d\tau' e^{(-i\nu' + i\omega)\tau'}
   \\ =
   \frac{1}{\beta^2} \sum_{n=-\infty}^\infty
   G(i\omega)
   \cdot \beta \delta_{\nu, \omega}
   \cdot \beta \delta_{\nu', \omega}
   \\ =
   \delta_{\nu, \nu'} G(i\nu)

Thus there is no scale factor relating the one and two frequency single particle Green's function

.. math::

   G(i\nu, i\nu') = \delta_{\nu, \nu'} G(i\nu)
