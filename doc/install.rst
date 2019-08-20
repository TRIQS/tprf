.. highlight:: bash

.. _install:

Packaged Versions of TPRF
=========================

A Docker image including the latest version of tprf is available `here <https://hub.docker.com/r/flatironinstitute/triqs>`_. For more information, please see the page on :ref:`TRIQS Docker <triqslibs:triqs_docker>`.

We further provide a Debian package for the Ubuntu LTS Versions 16.04 (xenial) and 18.04 (bionic), which can be installed by following the steps outlined :ref:`here <triqslibs:triqs_debian>`, and the subsequent command::

        sudo apt-get install -y triqs_tprf


Compiling TPRF from source
==========================

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

     $ source path_to_triqs/share/triqsvarsh.sh

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

    cmake ../tprf.src -DOPTION1=value1 -DOPTION2=value2 ... ../tprf.src

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
| Enable extra OpenMP threading                                   | -DCMAKE_CXX_FLAGS="-fopenmp -DTPRF_OMP"       |
+-----------------------------------------------------------------+-----------------------------------------------+
