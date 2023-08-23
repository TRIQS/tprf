(changelog)=

# Changelog

## Version 3.2.0

DESCRIPTION HERE, CLEANUP MESSAGES BELOW

We thank all contributors: Olivier Gingras, Alexander Hampel, Malte Rösner, Dylan Simon, Stefan Kaeser, Hugo U. R. Strand, Yann in 't Veld, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* Real-freq template of lindhard_chi00
* Fixed OMP thread safety lindhard_chi00
* Fixed thread safety of split_into_dynamic_wk_and_constant_k
* Added Eliashberg timings test script
* Remove use of deprecated TBLattice API Units
* No longer use periodization matrix but instead dims for mesh construction
* Port tprf to reworked triqs Green function meshes and interpolation
* Moved fermi and bose functions to lattice_utility.cpp
* simplify obtaining specific matsubara_freq in chi0_nr_from_gr_PH_at_specific_w
* Simplify tmesh boundary treatment
* Use const& more consistently
* Fix incorrect argument to gf operator[] in hubbard_atom.cpp
* Run TRIQS port_to_triqs3 script
* Use n_orb over orb_names in constructing hamiltonians
* Fix unused variable warnings
* Iterate on PR#32: Various smaller cleanups and improvements
* Run git clang-format on PR#32 commits
* Fixed documentation of split_into_dynamic_wk_and_constant_k function
* Fixed documentation in lindhard_chi00
* Updated documentation
* Added authorship to the headers
* Fix: Do not use floor division with complex floating point types
* Replace remaining occurances of gf_mesh + format
* Moved the duplicate function split_into_dynamic_and_static function out of gw.cpp and eliashberg.cpp into a seperate file
* Moved the dynamical_screened_interaction functions to a separate file
* Updated the lattice_desc file
* Made the freq-mesh argument notation consistent over all functions
* Put the real and imaginary freq. implementations of lindhard_chi00 under one overloaded function; adjusted the tests accordingly
* Made a test for the dynamical_screened_interaction_w function
* Updated documentation of lindhard_chi00_wk function; added chi00_square_lattice_fk to the tests to be run
* Adjusted lindhard_chi00_fk to match convention of band indices of lindhard_chi00_wk
* Fixed test for GW self-energy in real freq.
* Implementing tests for real freq. lindhard and gw sigma
* Implementation sigma_k_g0w0
* Implemented a GW routine for static interactions: needs testing
* Added documentation the the lindhard_chi00 functions
* Fixed the documentation of Matsubara frequency Green's function function
* OMP and MPI parallelized the eliashberg and lindhard_chi00 functions
* Add lattice_dyson_g_fk
* Add gw_sigma_fk_g0w0_spectral
* Add lattice_dyson_g0_fk
* Parallelize eliashberg_product
* Add orbital position to wannier90 reader
* Add Lindhard Chi0 in real frequencies
* move python vasp crpa parsers to tprf
* Update 3.1 easybuild script with sha256 of release tarball

### gf
* replace interp e_k(k) with indexing e_k[k]
* lattice_dyson_g_X possible OpenMP race condition fix
* Implemented a test for the lattice_dyson_g functions
* Made a template function for the lattice_dyson_g functions
* Initialized the Gfs to zero in the lattice_dyson_g functions
* Removed the enable extra OpenMP threading CMake option
* Implemented a few tests for the lattice_dyson_g_ functions; updated the gw test slightly
* Slight change of Gf initialization to 0
* fixed that lattice_dyson_g_ functions were not thread safe

### lattice_utils
* k_space_path -> triqs.lattice.utils.k_space_path

### util
* k_space_path rel and abs coord opt
* Templated the add_static_and_dynamic function

### doc
* k_space_path docstring update
* dbse tutorial wip
* sro tutorial add figs
* sro tutorial first complete draft
* sro dft w90 sc dmft
* sro tutorial viz scripts
* sro tutorial scaffold
* add new contributors
* theory add matrix rpa index warning
* theory chage cross and parallel chi to horiz and vert
* remove pyed deps in docs
* Corrected typos. Restructured some of the theory in the doc. Corrected error in the generation of the doc.
* Minimal changes to the ipynb for easier review.
* Updated tutorials to account for new TRIQS API and corrected typos.
* PH channel G diss typo
* Update Changelog for 3.1.1
* Add doc/_templates to sphinx exclude_patters, remove redundant doc/conf.py
* Provide more details on eliashberg, BSE and chi0 changes in changelog

### eli
* call tmesh.dlr_it().convolve_init() outside omp parallell region
* cppdlr convolve thread safety issue
* g delta g in DLR coeffspace
* Removed micro-timings
* Inconsistent handling of dlr and linear in preprocess_gamma_for_fft
* templated the g_delta_g_product function
* fix g delta g product using gf conjugate relation
* benchmark timings of g_delta_g
* dlr vs linear timing benchmark

### fitdlr
* fix bound and lstsq initial gues

### dlrfit
* stochastic tau fit routines

### test
* eli dlr solver tune cf decimals
* init mpi before tprf imports
* gw single kpts two band fix plot import
* fix gw test to gfv2 api
* lower eli dlr test runtimes
* Fix make_gf_dlr import in dlr_eliashberg_solver.py
* harden hf rpa test

### latt
* lattice_dyson_g_w do not store g_wk
* gf lattice_dyson_g_f doc
* gf real freq templated methods
* chi0 bubble in retime and g0 les gtr
* lindhard_chi00 real freq with fermi-dirac distrib not erfc
* enable real freq RPA

### dlr
* Added DLR to lindhard_chi00 function
* Added DLR to the remaining dynamical_screened_interaction_w functions
* rename dlr_coeff -> dlr
* Templated Eliashberg functions where possible
* Templated dlr_on_imtime and dlr_on_imfreq
* Eliashberg solver is working with DLR
* Added dlr_on_imtime and updated eliashberg test
* eliashberg_product_fft_constant is working
* Incorrect meshes in eliashberg_g_delta_g_product
* Eliashberg eliashberg_product_fft working
* Added dlr_on_imfreq for rank-4 gfs
* Eliashberg g_delta_g_product working
* fourier parallell fix
* gw hubbard dimer test passing
* gw routines working in gf.py test
* chi lattice bubble and fts
* gf fourier routines
* lattice dyson routines

### utils
* k_space_path accept list k input

### gw
* Updated documentation
* Added a single kpoint test for multiband system
* compute susceptibility
* fix dlr analytical continuation
* fix rho_k calc
* N_fix hdf5 tweak
* dlr solver working
* fix PI const
* Split the dynamic and static parts of the g0w_sigma function
* Renamed GW method comparison test
* small changes to GW compare methods test
* Added imaginary part in interaction of GW method comparison
* updated the documentation for index order fix
* Added a test which compares full GW and GW via spectral representation implementations
* Updated index order in real-freq GW
* save non int mu0
* extended single k-point gw test to real-freq
* fixed subtract V_k to get W_dyn_wk
* Added a single-kpoint gw test
* subtract V_k to get W_dyn_wk in gw_solver
* Hartree sigma V(q=0)rho(r=0)
* fix density calculations
* timers on post proc methods
* dynamical interaction fix index bug
* fix ed ref for hubbard dimer test
* parallelization optimization
* mpi print fixes
* add ase timer for benchmarking
* gw_dynamic_sigma renaming
* ret rf result
* hf real space calc and GWSolv h5 storage
* hubbard dimer ed reference result
* cleanup hubbard dimer tests
* move gw hubbard dimer to test
* spinless hubbard dimer test
* sc vs sic testing
* test cf ed, G0W0 and scGW
* GWSolver with tests
* hubbard dimer analytic ref test working
* analytic cf wip
* python doc for hartree_sigma fock_sigma
* hartree + fock sigma reprod hf_solver
* hf sigma tests wip
* Cleanup of g0w_sigma function
* Slightly shortened the g0w implementation
* Fixed inconsistent sign convention between static and dynamic parts
* Made a generic implementation for all the dynamical_screened interaction functions
* (Hopefully) fixed a bug in the gw static self-energy calculation which caused results that were sometimes wrong
* Added overloads to the dynamical_screened_interaction_W_from_generalized_susceptibility function
* a bit of cleanup
* fixed some bugs in the tests; added lindhard_chi00 to triqs_tprf.gw
* Added some documentation to the self-energy functions
* MPI parallelized the g0w_sigma function
* Removed the substraction of the static interaction in the dynamical_screened_interaction_w function; Updated the tests to account for this
* Updated the dynamical_screened_interaction_w tests
* Implemented tests for the GW self-energy functions
* overloaded functions for calculating the GW self-energy
* implemented calculation of GW self-energy with static interaction only; renamed some functions
* Added overloads to the dynamical_screened_interaction_W function for real-freq and dynamic interactions, and added some tests for these functions
* nda fixes and bugfix

### clang-tidy
* Synchronize config file with app4triqs

### numpy
* np.int -> np.int64

### gfv2
* port dlr parts to triqs/gfv2

### eliashberg
* Cleaned up the Eliashberg solver a bit

### dbse
* python solver
* c++ calculator

### g0w
* Renamed g0w_dyn_sigma
* Single kpoint for combined stat and dyn
* Implemented single kpoint calculation for g0w_dyn_sigma
* Single k-point calculations

### bse
* remove broken kmesh cmp
* attach_tri_vert

### sanitizer
* error fixes

### bubble
* verbose flag

### local
* chi0 tau bubble from g tau

### fourier
* r<->k fft for real freq and real time gfs

### cmake
* Bump Version-number to 3.1.1

### py
* one more numpy depr fix
* fix numpy depr warnings (errors with modern numpy)

### lat
* PH vector product


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
