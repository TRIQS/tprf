..
   Generated automatically by cpp2rst

.. highlight:: c


.. _single_particle_greens_function:

single_particle_greens_function
===============================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_iw_t single_particle_greens_function (int nw, double beta, double U)

Single-particle Green's function of the Hubbard atom at half-filling


Parameters
----------

 * **nw**: number of Matsubara frequencies

 * **beta**: inverse temperature

 * **U**: Hubbard U interaction parmeter



Return value
------------

g single-particle Green's function of the Hubbard atom :math:`G(i\omega_n)`

Documentation
-------------


     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
     please cite the paper if you use this function!

     .. math::
         G(i\omega_n) = \frac{1}{i\omega_n - \frac{U^2}{4 i\omega_n}}