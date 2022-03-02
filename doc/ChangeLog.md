(changelog)=

# Changelog

## Version 3.1.0

TPRF version 3.1.0 is a compatibility release for TRIQS version 3.1.0 also containing some new functionality.

Eliashberg
----------
* Functionality for solving the linearized Eliasberg equation

Bethe Salpeter Equation (BSE)
-----------------------------
* Functionality to solve the BSE at one (finite) bosonic frequency at a time

Contributors: Stefan Käser, Hugo U.R. Strand

Maintenance 
-----------
* Change np.complex -> complex
* Compiler warning fixes
* Thread race-condition bugfix

TRIQS compatibility
-------------------
* Compatibility updates to adhere to changes in TRIQS from v3.0.0 to v3.1.0
* Merges from app4triqs to adhere to changes in the TRIQS application framework from v3.0.0 to v3.1.0

Contributors: Nils Wentzell, Alexander Hampel, Dylan Simon, Hugo U.R. Strand 

## Version 3.0.0

tprf version 3.0.0 is a compatibility
release for TRIQS version 3.0.0 that
* introduces compatibility with Python 3 (Python 2 no longer supported)
* adds a cmake-based dependency management
* fixes various smaller application issues

We provide a more detailed description of the changes below.

General
-------
* Protect various solver logos for non-utf encoding of sys.stdout
* Rename pytriqs->triqs
* Run port_to_triqs3 script
* Merge app4triqs python3 changes
* Fix bug in wannier90 hr reader

cmake
-----
* Use find_package(OpenMP ...) and define openmp INTERFACE target
* Link cpp2py module also against triqs_py library

doc
---
* Fix sidebar version number
* Add a section on Anaconda to the install page

py3
---
* Use 2to3 to port python files and notebooks
* Fix floor division in various places
* Make sure to use spaces for indentation and no tabs

python
------
* Make sure to import pytriqs.utility.mpi in lattice module init

Contributors: Philipp Dumitrescu, Dylan Simon, Nils Wentzell, Manuel Zingl

## Version 2.2.0

TPRF version 2.2.0 is a compatibility release
for TRIQS version 2.2.0. It provides improvements to
the documentation and fixes various smaller issues.

We provide a more detailed description of the changes below.

doc
---
* Include debian package on installation page
* Correct triqs links in conf.py.in
* Fixes of various spelling error
* Remove generated documentation, only regenerate on doc build

fourier
-------
* Fix sanitizer positives related to triqs unstable view adjustments

General
-------
* FIX in bse.py, We cannot generally guarantee that Idx(0,1,2) is a valid index of the mesh
* FIX Do not use bracket operator of gf when domain_pt of cluster mesh is passed
* Instead of gf<..>::zero_t use ::target_t::value_t
* Use zeros(g.target_shape()) instead of g.get_zero()
* Changes to restore triqs/2.2.0 compatibility

hf
--
* Explicitly check compatibility of H_int and gf_struct


## Version 2.1.1

* We now provide a debian packaged version of tprf which is also part of the [triqs docker image](https://hub.docker.com/r/flatironinstitute/triqs)
* Updated documentation with debian package install instructions
* Minor fixes in the Documentation pages (spelling, corrected links)
* Added a check for the compatibility between H_int and gf_struct
