.. _dbse_hubbard_zeeman:

Susceptibility of the Hubbard-Zeeman model
==========================================

In this tutorial we will compute the dynamic, momentum-resolved susceptibility :math:`\chi_{abcd}(\omega,\mathbf{q})` of the Hubbard-Zeeman model, along the lines of `Communications Physics 6, 289 (2023) <https://doi.org/10.1038/s42005-023-01411-w>`_. Unlike previous tutorials, in this tutorial we will compute and look at the full tensor structure of the susceptibility. Both the BSE and DBSE will be used. 

The model has a hopping parameter :math: `t`, a Hubbard interaction :math: `U`, and a magnetic field :math: `B`. The Zeeman field lifts the degeneracy between spin up and spin down and leads to a richer spin structure of the susceptibility including non-trivial off-diagonal elements. TRIQS/cthyb is used as an impurity solver for the DMFT self-consistency loop, and then the vertices are calculated using W2Dynamics. The overall structure of the scripts is similar to previous tutorials. (insert scripts here)


Vertices
--------

(plotting scripts for vertices here)


Bethe-Salpeter Equation
-----------------------

(short discussion of BSE script)

Convergence of Susceptibility
-----------------------------

(Show convergence of different components here)

Larmor precession
-----------------

In the presence of a Zeeman field :math: `B`, the total spin of the system will undergo Larmor precession. This is visible in the in-plane (x,y) components of the spin susceptibility. The frequency of the oscillation is fixed by :math: `B` and is thus known analytically. This is a good check on the numerical results. 

(plot scripts)