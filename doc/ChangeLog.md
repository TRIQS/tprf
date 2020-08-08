Version 3.0.0
=============

tprf version 3.0.0 is a compatibility
release for TRIQS version 3.0.0 that
* introduces compatibility with Python 3 (Python 2 no longer supported)
* adds a cmake-based dependency management
* fixes several application issues


Version 2.2.0
=============

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


Version 2.1.1
-------------

* We now provide a debian packaged version of tprf which is also part of the [triqs docker image](https://hub.docker.com/r/flatironinstitute/triqs)
* Updated documentation with debian package install instructions
* Minor fixes in the Documentation pages (spelling, corrected links)
* Added a check for the compatibility between H_int and gf_struct
