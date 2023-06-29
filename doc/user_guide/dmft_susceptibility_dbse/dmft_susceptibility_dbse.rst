.. _dmft_susceptibility_dbse:

Spin susceptibility in Sr2RuO4
==============================

In this tutorial we will compute the static magnetic susceptibilitly :math:`\chi_{S_z S_z}(\mathbf{q})` of the correlated Hund's metal Sr2RuO4 withn dynamical mean field theory (DMFT), reproducing the results of https://doi.org/10.1103/PhysRevB.100.125120. We will use the reformulation of the (DMFT) lattice susceptibility using dual propagators, for details see https://arxiv.org/abs/2306.05157.

The calculation is based on a Wannier model for the three bands crossing the Fermi level in Sr2RuO4. These bands have Ru-4d t2g symmetry and a Wannier interpolation with Wannier90 converges in just a few iterations, giving the band structure

.. image:: figure_sro_band_structure.svg
    :align: center

see :download:`tight_binding_model.py <tight_binding_model.py>`.

The Wannier Hamiltonian is combined with a local Kanamori interaction with Hubbard :math:`U=2.4` eV and Hund's :math:`J=0.4` eV and the self-consistent DMFT solution is determined using TRIQS/cthyb as impurity solver. The scripts for the DMFT solution are :download:`common.py <common.py>` and :download:`calc_sc_dmft.py <calc_sc_dmft.py>`.


Wannier bandstructure
---------------------

.. literalinclude:: tight_binding_model.py
   :lines: 23-

.. literalinclude:: plot_band_structure.py
   :lines: 23-

.. image:: figure_sro_band_structure.svg
    :align: center

Impurity correlators from W2Dynamics
------------------------------------

.. literalinclude:: calc_g2.py
   :lines: 23-

.. literalinclude:: calc_tri.py
   :lines: 23-

.. literalinclude:: calc_chi.py
   :lines: 23-

Dual Bethe-Salpeter equation
----------------------------

.. literalinclude:: calc_dbse.py
   :lines: 23-
   
.. image:: figure_sro_chi_bandpath.svg
   :align: center
   :width: 400

Plot script :download:`plot_dbse.py <plot_dbse.py>`

	   
Old
-------------------------------------------------------

In this guide we will compute the uniform magnetic susceptibility :math:`\chi = \chi(\mathbf{Q} = \mathbf{0})` of the single band Hubbard model on the square lattice with nearest neighbour hopping using dynamical mean-field theory (DMFT).

We will do this in two very different ways from

1. self consistent DMFT calculations in applied fields, and from
2. the lattice Bethe-Salpeter Equation using the DMFT local vertex.

Since DMFT is thermodynamically consistent *[Hafermann et al., PRB 90, 235105 (2014)]* these two approaches gives the same susceptibility within error-bars.

Lattice susceptibility from the Bethe-Salpeter Equation
-------------------------------------------------------

Instead of multiple calculations in applied field the susceptibility :math:`\chi` can be obtained from a direct calculation of the two-particle susceptbility using the DMFT local vertex and the lattice Bethe-Salpeter equation.

To this end one has to perform the steps

1. compute the DMFT impurity two-particle Green's function :math:`G^{(2)}`,
2. compute the DMFT impurity two-particle vertex :math:`\Gamma`, and
3. solve the lattice Beht-Salpeter Equation for the lattice susceptibility :math:`\chi(\mathbf{Q})`.

DMFT local vertex
^^^^^^^^^^^^^^^^^

To obtain the DMFT impurity magnetic vertex :math:`\Gamma_m` one first computes the DMFT impurity two-particle Green's function in the particle-hole channel :math:`G^{(2,PH)}_{abcd}(i\omega, i\nu, i\nu')`, where :math:`a,b,c,d` are spin-orbital indices. The two-particle Green's function can be sampled using `triqs_cthyb` and to compute the static susceptibility it is sufficient to only keep the zero Bosonic frequency :math:`\omega = 0`.

For the single band model the magnetic susceptbility :math:`\chi_m` is directly related to :math:`G^{(2,PH)}` as

.. math::
   \chi_m(i\omega, i\nu, i\nu') =
   G^{(2,PH)}_{\uparrow\uparrow\uparrow\uparrow}(i\omega, i\nu, i\nu')
   -
   G^{(2,PH)}_{\uparrow\uparrow\downarrow\downarrow}(i\omega, i\nu, i\nu')
   \, ,

and the magnetic bubble susceptibility :math:`\chi^{(0)}_m` is given by :math:`\chi^{(0)}_m(i\omega. i\nu, i\nu') = - \beta \delta_{\nu, \nu'} G(i\nu) G(i\omega + i\nu)`. The two susceptibilities are related by the impurity Bethe-Salpeter Equation giving the corresponding magetic vertex function :math:`\Gamma_m` as

.. math::
   \Gamma_m = [\chi^{(0)}_m]^{-1} - [\chi_m]^{-1}

Starting from the self consistent DMFT solution at zero-field we run `triqs_cthyb` to sample the two-particle Green's function and then compute :math:`\chi_m`, :math:`\chi_m^{(0)}`, and :math:`\Gamma_m` see below
   
.. literalinclude:: calc_g2.py
   :lines: 23-

The resulting response functions are plotted below	   

.. image:: figure_g2.svg
    :align: center

The visualization script is available here: :download:`plot_g2.py <plot_g2.py>`.

Lattice Bethe-Salpeter Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equipped with the DMFT local vertex :math:`\Gamma_m` it is possible to compute the DMFT lattice susceptibility :math:`\chi(\mathbf{Q})` from the lattice Bethe-Salpeter Equation (BSE)

.. math::
   \chi(\mathbf{Q}) = \chi_0(\mathbf{Q}) - \chi_0(\mathbf{Q}) \Gamma \chi(\mathbf{Q})
   \, .

TPRF comes with an OpenMP amd MPI parallelized BSE solver `triqs_tprf.bse.solve_lattice_bse`. However, the calculation is done with a fixed number of frequencies :math:`n_\nu` in the fermionic frequencies :math:`\nu` and :math:`\nu'`, and the solution converges only linearly with the size of the frequency window. Therefore we solve the BSE for a range of window sizes to enable extrapolation :math:`N_\nu \rightarrow \infty`.
   
.. literalinclude:: calc_bse.py
   :lines: 23-

The resuls along the high symmetry path of the Brillouin zone is shown below for fixed :math:`N_\nu` (left panel) and the extrapolation for the :math:`\Gamma`-point is also shown (right panel).
	   
.. image:: figure_bse.svg
    :align: center

The visulaization script is available here: :download:`plot_bse.py <plot_bse.py>`.

The result for the homogeneous magnetic susceptbilitiy :math:`\chi(\mathbf{0})` from the BSE is

.. math::
   \chi_{\textrm{BSE}} = \lim_{N_\nu \rightarrow \infty} \chi(\mathbf{0}) \approx 0.3472

in quantitative agreement with the applied field value.
   
Summary
-------

Now we can compare the two results for the homogeneous static magnetic susceptibility, from

1. the applied field calculation, and
2. the BSE calculation.

The results are in quantitative agreement

.. math::
   \chi_{\textrm{Field}} \approx 0.3479
   \, ,

.. math::
   \chi_{\textrm{BSE}} \approx 0.3472
   \, ,

and the accuracy is limited by the stochastic Monte Carlo noise in the vertex. This can be improved further by increasing the number of samples in the `triqs_cthyb` two-particle Green's function calculation.

Note, that the BSE approach is much more general than the applied field approach. The BSE calculation gives the whole momentum dependent lattice susceptibility :math:`\chi(\mathbf{Q})` and provides dynamical information when using finite Bosonic freqiencies :math:`|\omega| > 0`.
