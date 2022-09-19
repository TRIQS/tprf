(changelog)=

# Changelog

## Version 3.1.1

TPRF version 3.1.1 is a patch release that fixes an issue with recent numpy versions.

We thank all contributors: Stefan, Hugo U. R. Strand, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* Update 3.1 easybuild script with sha256 of release tarball

### py
* Fix numpy depr warnings (errors with modern numpy)

### doc
* Provide more details on eliashberg, BSE and chi0 changes in changelog


## Version 3.1.0

TPRF version 3.1.0 is a compatibility release for TRIQS version 3.1.0 also containing some new functionality.

Eliashberg
----------
* Added functionality to solve the multi-orbital linearized Eliashberg equation for the most dominant gap functions:
    * Only irreducible particle-particle vertices dependent on one bosonic frequency and one momentum are supported.
    * Matrix product inside Eliashberg equation is implemented in two ways:
        * Explicit loops
        * Fourier transformations taking advantage of the convolution theorem
     Depending on the used meshs either implementation can be more efficient.
    * The Fourier implementation is parallelized over threads.
    * Functionality to solve for the most dominant gap function using the power method or solving for the `k` most dominant ones using the implicitly restarted Arnoldi method.
    * Functionality to only solve for gap functions with specific symmetries (frequency, momentum or orbital).
    * Added constructors for the singlet and triplet irreducible vertices in the random phase approximation.

Bethe-Salpeter equation (BSE)
-----------------------------
* Added functionality to solve the BSE at specific bosonic frequencies individually to avoid memory limitations.

Chi0
----
* Added memory optimized version of chi0 construction for smaller frequency meshes.
* Added functionality to construct chi0 for specific frequency.

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
