.. highlight:: bash

.. _install:

Packaged Versions of TPRF
=========================

.. _ubuntu_debian:
Ubuntu Debian packages
----------------------

We provide a Debian package for the Ubuntu LTS Versions 16.04 (xenial) and 18.04 (bionic), see :ref:`Triqs Debian package install instructions <triqslibs:ubuntu_debian>`. After installing Triqs the command::

  sudo apt-get install -y triqs_tprf

can be used to install the TPRF package.
    
.. _anaconda:
Anaconda (experimental)
-----------------------

We provide Linux and OSX packages for the `Anaconda <https://www.anaconda.com/>`_ distribution. The packages are provided through the `conda-forge <https://conda-forge.org/>`_ repositories. After `installing conda <https://docs.conda.io/en/latest/miniconda.html>`_ you can install TPRF with::

     conda install -c conda-forge triqs_tprf

See also `github.com/conda-forge/triqs_tprf-feedstock <https://github.com/conda-forge/triqs_tprf-feedstock/>`_.

.. _docker:
Docker
------

A Docker image including the latest version of TPRF is available `here <https://hub.docker.com/r/flatironinstitute/triqs>`_. For more information, please see the page on :ref:`TRIQS Docker <triqslibs:triqs_docker>`.


Compiling TPRF from source
==========================

.. note:: To guarantee reproducibility in scientific calculations we strongly recommend the use of a stable `release <https://github.com/TRIQS/triqs/releases>`_ of both TRIQS and its applications.

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/tprf`` repository from GitHub::

     $ git clone https://github.com/TRIQS/tprf tprf.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir tprf.build && cd tprf.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqs/triqsvars.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../tprf.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Version compatibility
---------------------

Keep in mind that the version of ``tprf`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd tprf.src && git tag

Checkout the version that you want::

     $ git checkout 2.1.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``tprf`` can be configured using CMake-options::

    cmake ../tprf.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_tprf           |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
