.. _dbse_hubbard_zeeman:

Susceptibility of the Hubbard-Zeeman model
==========================================

In this tutorial we will compute the dynamic, momentum-resolved susceptibility :math:`\chi_{abcd}(\omega,\mathbf{q})` of the Hubbard-Zeeman model, along the lines of `Communications Physics 6, 289 (2023) <https://doi.org/10.1038/s42005-023-01411-w>`_. Unlike previous tutorials, in this tutorial we will compute and look at the full tensor structure of the susceptibility. Both the BSE and DBSE will be used. 

The model has a hopping parameter :math: `t`, a Hubbard interaction :math: `U`, and a Zeeman magnetic field :math: `B`. The latter appears in the Hamiltonian as :math:`-B (\hat{n}_\uparrow-\hat{n}_\downarrow)` The Zeeman field lifts the degeneracy between spin up and spin down and leads to a richer spin structure of the susceptibility including non-trivial off-diagonal elements. TRIQS/cthyb is used as an impurity solver for the DMFT self-consistency loop, and then the vertices are calculated using W2Dynamics. The overall structure of the scripts is similar to previous tutorials: :download:`common.py <common.py>` is a module with general functionality, :download:`calc_sc_dmft.py <calc_sc_dmft.py>` runs the DMFT self-consistency loop and then :download:`calc_g2.py <calc_g2.py>`, :download:`calc_tri.py <calc_tri.py>`, :download:`calc_chi.py <calc_chi.py>` measure the three required correlation functions. Finally, :download:`calc_bse.py <calc_bse.py>` evaluates the Bethe-Salpeter equation with different fermionic cut-offs ``nwf``. 

Vertices
--------

:download:`plot_vertex.py <plot_vertex.py>`
(Include discussion of plot script and figure here)

Convergence of Susceptibility
-----------------------------

:download:`plot_components.py <plot_components.py>`

(Include discussion of plot script and figure here)

Larmor precession
-----------------

In the presence of a Zeeman field :math:`B`, the total spin of the system will undergo Larmor precession. This is visible in the in-plane (x,y) components of the spin susceptibility at :math:`\mathbf{q}=\Gamma`. The frequency of the oscillation is fixed by :math:`B` and is thus known analytically. This is a good check on the numerical results. 

:download:`plot_chiG.py <plot_chiG.py>`

(Include discussion of plot script and figure here)