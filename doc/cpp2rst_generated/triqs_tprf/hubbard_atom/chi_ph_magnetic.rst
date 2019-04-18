..
   Generated automatically by cpp2rst

.. highlight:: c


.. _chi_ph_magnetic:

chi_ph_magnetic
===============

**Synopsis**:

.. code-block:: c

    triqs_tprf::g2_iw_t chi_ph_magnetic (int nw, int nwf, double beta, double U)

Magnetic susceptibility of the Hubbard atom at half-filling :math:`\chi(\omega, \nu, \nu')`


Parameters
----------

 * **nw**: number of bosonic Matsubara frequencies

 * **nwf**: number of fermionic Matsubara frequencies

 * **beta**: inverse temperature

 * **U**: Hubbard U interaction parmeter



Return value
------------

chi magnetic susceptibility :math:`\chi(\omega, \nu, \nu')`

Documentation
-------------


     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
     please cite the paper if you use this function!