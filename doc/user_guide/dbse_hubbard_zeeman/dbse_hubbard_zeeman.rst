.. _dbse_hubbard_zeeman:

Susceptibility of the Hubbard-Zeeman model
==========================================

In this tutorial we will compute the dynamic, momentum-resolved susceptibility :math:`\chi_{abcd}(\omega,\mathbf{q})` of the Hubbard-Zeeman model, along the lines of `Communications Physics 6, 289 (2023) <https://doi.org/10.1038/s42005-023-01411-w>`_. Unlike previous tutorials, in this tutorial we will compute and look at the full tensor structure of the susceptibility. Both the BSE and DBSE will be used. 

The model has a hopping parameter :math: `t`, a Hubbard interaction :math: `U`, and a Zeeman magnetic field :math: `B`. The latter appears in the Hamiltonian as :math:`-B (\hat{n}_\uparrow-\hat{n}_\downarrow)` The Zeeman field lifts the degeneracy between spin up and spin down and leads to a richer spin structure of the susceptibility including non-trivial off-diagonal elements. TRIQS/cthyb is used as an impurity solver for the DMFT self-consistency loop, and then the vertices are calculated using W2Dynamics. The overall structure of the scripts is similar to previous tutorials: :download:`common.py <common.py>` is a module with general functionality, :download:`calc_sc_dmft.py <calc_sc_dmft.py>` runs the DMFT self-consistency loop and then :download:`calc_g2.py <calc_g2.py>`, :download:`calc_tri.py <calc_tri.py>`, :download:`calc_chi.py <calc_chi.py>` measure the three required correlation functions. Finally, :download:`calc_bse.py <calc_bse.py>` evaluates the Bethe-Salpeter equation with different fermionic cut-offs ``nwf``. 

Convergence of Susceptibility
-----------------------------

The following plot script shows the convergence with the fermionic frequency box size of the normal and dual Bethe-Salpeter equation, in color and black respectively. We focus on one particular component of the orbital tensor and one bosonic frequency :math: `\Omega=3`. The linear scaling with :math: `1/N_\nu` for the normal Bethe-Salpeter equation is clearly visible, which makes extrapolation necessary. On the other hand, for the dual Bethe-Salpeter equation, the calculations with finite frequency box size are very close to the extrapolated DBSE result. 

In this plot, the extrapolated values of the BSE and DBSE are not identical. The BSE and DBSE converge to the same result once (1) the Monte Carlo error bars on all impurity Green's functions are negligible (2) the DMFT self-consistency loop is converged perfectly (3) the fermionic box sizes are sufficiently large that extrapolation can be done. In practice, given a finite amount of computational resources for the Monte Carlo, the two results are expected to be close.

:download:`plot_components.py <plot_components.py>`

.. image:: figure_bse_w0_comp12.svg
   :align: center

Larmor precession
-----------------

In the presence of a Zeeman field :math:`B`, the total spin of the system will undergo Larmor precession. This is visible in the in-plane (x,y) components of the spin susceptibility at :math:`\mathbf{q}=\Gamma`. The frequency of the oscillation is fixed by :math:`B` and is thus known analytically. This is a good check on the numerical results. 

:download:`plot_chiG.py <plot_chiG.py>`

.. image:: figure_bse_w0_comp12.svg
   :align: center

From left to right, this plot shows the :math: `S^x S^x`, :math: `S^x S^y`, :math: `S^z S^z` and :math: `N N` components of :math: `\chi(q=\Gamma,\omega_m)` as a function of :math: `\omega_m`. The orange lines are the BSE result, which darker lines indicating a larger fermionic frequency box. The blue lines are the DBSE results for the same frequency box sizes, which show essentially perfect convergence. In the first two plots, the black lines are the analytically known exact result. For more information, see `Communications Physics 6, 289 (2023) <https://doi.org/10.1038/s42005-023-01411-w>`_.