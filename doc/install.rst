.. highlight:: bash

.. _install:

Compiling TPRF from source
==========================

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that Triqs is installed in the ``path_to_triqs`` directory.

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

Be careful that the version of the TRIQS library and of the application must be
compatible (more information on the :ref:`TRIQS website <triqslibs:versions>`).
In particular you should make sure that the Major and Minor Version number
of the application and TRIQS agree.
If you want to use a particular version of the application, go into the directory with the sources
and look at all available versions::

     $ cd tprf.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

Then follow the steps 2 to 4 described above to compile the code.

Custom CMake options
--------------------

The functionality of ``tprf`` can be tweaked using extra compile-time options passed to CMake::

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
