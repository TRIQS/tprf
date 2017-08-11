.. highlight:: bash

.. _install:

Installation
============


Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

Installation steps
------------------

#. Download the sources of the solver from github::

     $ git clone https://github.com/HugoStrand/tprf.git tprf.src

#. Create an empty build directory where you will compile the code::

     $ mkdir tprf.build && cd tprf.build

#. In the build directory call cmake specifying where the TRIQS library is installed::

     $ cmake -DTRIQS_PATH=path_to_triqs ../tprf.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

.. note:: Be careful with the cmake command above: set TRIQS_PATH, not CMAKE_INSTALL_PREFIX (this variable is only for the TRIQS library)!

Version compatibility
---------------------

Be careful that the version of the TRIQS library and of the solver must be
compatible (more information on the :ref:`TRIQS website <triqslibs:versions>`).
If you want to use a version of
the solver that is not the latest one, go into the directory with the sources
and look at all available versions::

     $ cd tprf.src && git tag

Checkout the version of the code that you want::

     $ git checkout 1.0.0

Then follow the steps 2 to 4 described above to compile the code.

Custom CMake options
--------------------

Functionality of ``tprf`` can be tweaked using extra compile-time options passed to CMake::

    cmake -DTRIQS_PATH=path_to_triqs -DOPTION1=value1 -DOPTION2=value2 ... ../cthyb.src

+---------------------------------------------------------------+-----------------------------------+
| Options                                                       | Syntax                            |
+===============================================================+===================================+
| Disable testing (not recommended)                             | -DTests=OFF                       |
+---------------------------------------------------------------+-----------------------------------+
| Build the documentation locally                               | -DBUILD_DOC=ON                    |
+---------------------------------------------------------------+-----------------------------------+
