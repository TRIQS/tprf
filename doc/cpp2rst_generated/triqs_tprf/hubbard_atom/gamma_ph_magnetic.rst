..
   Generated automatically by cpp2rst

.. highlight:: c


.. _gamma_ph_magnetic:

gamma_ph_magnetic
=================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g2_iw_t gamma_ph_magnetic (int nw, int nwf, double beta, double U)

Magnetic vertex function in the particle-hole channel of the Hubbard atom at half-filling :math:`\Gamma(\omega, \nu, \nu')`


Parameters
----------

 * **nw**: number of bosonic Matsubara frequencies

 * **nwf**: number of fermionic Matsubara frequencies

 * **beta**: inverse temperature

 * **U**: Hubbard U interaction parmeter



Return value
------------

gamma magnetic susceptibility :math:`\Gamma(\omega, \nu, \nu')`

Documentation
-------------


    Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
    please cite the paper if you use this function!