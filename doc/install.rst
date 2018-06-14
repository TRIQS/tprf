.. highlight:: bash

.. _install:

Installation
============


Prerequisite
-------------------

#. The TRIQS library (https://github.com/TRIQS/triqs) and its dependecies. See the install instructions in the TRIQS documentation (https://triqs.github.io/triqs/).

Installation steps
------------------

#. Clone the github repository::

     $ git clone https://github.com/HugoStrand/tprf.git tprf.src

#. Create an empty build directory where you will compile the code::

     $ mkdir tprf.build && cd tprf.build

#. In the build directory call cmake::

     $ cmake ../tprf.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make install
     $ make test

Version compatibility
---------------------

Be careful that the version of the TRIQS library and of the solver must be
compatible.

Custom CMake options
--------------------

The functionality of ``tprf`` can be tweaked using extra compile-time options passed to CMake::

    cmake -DOPTION1=value1 -DOPTION2=value2 ... ../cthyb.src

+---------------------------------------------------------------+-----------------------------------+
| Options                                                       | Syntax                            |
+===============================================================+===================================+
| Build the documentation locally                               | -DBUILD_DOC=ON                    |
+---------------------------------------------------------------+-----------------------------------+
