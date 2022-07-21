.. highlight:: python

.. _gf_functions: 

Lattice Green's functions
=========================

.. autofunction:: triqs_tprf.lattice.lattice_dyson_g0_wk
.. autofunction:: triqs_tprf.lattice.lattice_dyson_g0_fk
.. autofunction:: triqs_tprf.lattice.lattice_dyson_g_wk
.. autofunction:: triqs_tprf.lattice.lattice_dyson_g_fk
.. autofunction:: triqs_tprf.lattice.lattice_dyson_g_w
		  
Lindhard non-interacting generalized susceptibility
===================================================

.. autofunction:: triqs_tprf.lattice.lindhard_chi00

Random Phase Approximation
==========================

.. autofunction:: triqs_tprf.lattice.solve_rpa_PH
.. autofunction:: triqs_tprf.rpa_tensor.kanamori_quartic_tensor

Impurity susceptibility and Bethe-Salpeter Equation
===================================================

.. autofunction:: triqs_tprf.chi_from_gg2.chi0_from_gg2_PH
.. autofunction:: triqs_tprf.chi_from_gg2.chi0_from_gg2_PP
.. autofunction:: triqs_tprf.chi_from_gg2.chi_from_gg2_PH
.. autofunction:: triqs_tprf.chi_from_gg2.chi_from_gg2_PP

.. autofunction:: triqs_tprf.bse.solve_local_bse
   
Lattice Bethe-Salpeter Equation
===============================

.. autofunction:: triqs_tprf.bse.solve_lattice_bse
.. autofunction:: triqs_tprf.bse.solve_lattice_bse_at_specific_w
.. autofunction:: triqs_tprf.bse.get_chi0_wnk
.. autofunction:: triqs_tprf.bse.get_chi0_nk_at_specific_w
		  
GW approximation
================

.. autofunction:: triqs_tprf.lattice.dynamical_screened_interaction_W
.. autofunction:: triqs_tprf.lattice.dynamical_screened_interaction_W_from_generalized_susceptibility
.. autofunction:: triqs_tprf.gw.bubble_PI_wk
.. autofunction:: triqs_tprf.gw.gw_sigma
.. autofunction:: triqs_tprf.gw.g0w_sigma

Linearized Eliashberg equation 
==============================

.. autofunction:: triqs_tprf.eliashberg.solve_eliashberg
.. autofunction:: triqs_tprf.eliashberg.preprocess_gamma_for_fft
.. autofunction:: triqs_tprf.eliashberg.semi_random_initial_delta
.. autofunction:: triqs_tprf.eliashberg.power_method_LR
.. autofunction:: triqs_tprf.eliashberg.implicitly_restarted_arnoldi_method
.. autofunction:: triqs_tprf.eliashberg.construct_gamma_singlet_rpa
.. autofunction:: triqs_tprf.eliashberg.construct_gamma_triplet_rpa

Hubbard atom analytic response functions
========================================

.. autofunction:: triqs_tprf.analytic_hubbard_atom.analytic_hubbard_atom

Two-particle response function linear-algebra operations
========================================================

.. autofunction:: triqs_tprf.linalg.inverse_PH
.. autofunction:: triqs_tprf.linalg.inverse_PP
.. autofunction:: triqs_tprf.linalg.inverse_PH_bar

.. autofunction:: triqs_tprf.linalg.product_PH
.. autofunction:: triqs_tprf.linalg.product_PP
.. autofunction:: triqs_tprf.linalg.product_PH_bar

.. autofunction:: triqs_tprf.linalg.identity_PH
.. autofunction:: triqs_tprf.linalg.identity_PP
.. autofunction:: triqs_tprf.linalg.identity_PH_bar
		  
Wannier90 tight binding parsers
===============================

.. autofunction:: triqs_tprf.wannier90.parse_hopping_from_wannier90_hr_dat
.. autofunction:: triqs_tprf.wannier90.parse_lattice_vectors_from_wannier90_wout
.. autofunction:: triqs_tprf.wannier90.parse_reciprocal_lattice_vectors_from_wannier90_wout
.. autofunction:: triqs_tprf.wannier90.parse_band_structure_from_wannier90_band_dat

Tight binding lattice model
===========================

.. autoclass:: triqs_tprf.tight_binding.TBLattice
   :members:

.. autofunction:: triqs_tprf.tight_binding.create_square_lattice

.. autoclass:: triqs_tprf.super_lattice.TBSuperLattice
   :members:
      
Hartree-Fock and Hartree solvers
================================

.. autoclass:: triqs_tprf.hf_solver.HartreeFockSolver
   :members:

.. autoclass:: triqs_tprf.hf_response.HartreeFockResponse
   :members:

.. autoclass:: triqs_tprf.hf_solver.HartreeSolver
   :members:
      
.. autoclass:: triqs_tprf.hf_response.HartreeResponse
   :members:

Parameter collections
=====================

.. autoclass:: triqs_tprf.ParameterCollection.ParameterCollection
   :members:
.. autoclass:: triqs_tprf.ParameterCollection.ParameterCollections
   :members:
.. autofunction:: triqs_tprf.ParameterCollection.parameter_scan
