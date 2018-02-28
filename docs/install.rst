Installing hophop
=================

To compile and run *hophop*, you need the following requirements:

Requirements
------------
The following libraries and software is required to compile successfully:

* C compiler
    Any compiler fulfilling the C99 standard is suitable.
    *hophop* was tested with `Gnu Compiler Collection <https://gcc.gnu.org/>`_
    and `LLVM CLang <https://clang.llvm.org/>`_ (with no significant speed
    differences).
* `CMake <https://cmake.org/>`_ > 3.0
    We use CMake as a build system, so you need to have it installed.
* `Gnu Scientific Library (GSL) <https://www.gnu.org/software/gsl/>`_
    The gnu scientific library provides many mathematical algorithms,
    functions, constants and so on. It is heavily used within the program,
    mainly for pseudo random number generation and probability
    distributions.
* `OpenMP <http://www.openmp.org/>`_
    Required for shared-memory parallelization. This is usually shipped with 
    the compiler.
* `Lis <http://www.ssisc.org/lis/>`_ >= 1.4.43 (optional, but recommended!)
    Lis (Library of Iterative Solvers for linear systems) is used as a solver
    for the balance equations method. When it is not found, the `mgmres` solver
    is used, the source code of which is shipped with `hophop`. However, we
    recommend using Lis where possible.

.. note:: You may find some of the requirements in the repositories of your Linux distribution, at least the compiler,
          CMake, and OpenMP. On Debian or Ubuntu Linux, for example, you can simply run the following command
          to download and install some the requirements: ::

              $ apt-get install build-essential cmake

Downloading the code
--------------------

Please either clone the `master` branch from `hophop's Github repository <https://github.com/janoliver/hophop>`_
or download one of the stable releases from `hophop's Release page <https://github.com/janoliver/hophop/releases>`_.

Building
--------

With all the requirements in standard (i.e., discoverable by CMake) paths,
you may be lucky and the following works instantly: ::

    $ tar xzf hophop-2.4.tar.gz
    $ mkdir build_hophop
    $ cd build_hophop
    $ cmake ../hophop-2.4
    $ make
    $ make install


.. Tip:: You can change the `install` location with the :code:`CMAKE_INSTALL_PREFIX` 
         command line variable: ::

             $ cmake ../hophop-2.4 -DCMAKE_INSTALL_PREFIX=/usr/local

When CMake can't figure out the locations of `Lis` and `GSL`, you can specify the
following variables to help searching:

* :code:`LIS_ROOT_DIR=/path/to/lis/location`
* :code:`GSL_ROOT_DIR=/path/to/gsl/location`

The locations must contain an :code:`include/` and :code:`lib` (or :code:`lib64`)
folder, where the headers and libraries are located. Example: ::

    $ cmake ../hophop-2.4 -DLIS_ROOT_DIR=/opt/lis/ -DGSL_ROOT_DIR=/opt/gsl
    $ make
    $ make install

.. Tip:: You can save custom library locations in the binary's `rpath`, so they
         are found without requiring :code:`LD_LIBRARY_PATH` to be set. ::

             $ cmake ../hophop-2.4 -DCMAKE_EXE_LINKER_FLAGS='-Wl,-rpath,/opt/lis/lib:/opt/gsl/lib'

Running hophop
--------------

When everything is built correctly, you can try running *hophop* by simply typing ::

    $ /path/to/install/location/hophop

It will use some default parameters to run. 

.. Tip:: When you get some :code:`library not found` errors, set the 
         :code:`LD_LIBRARY_PATH` variable to the location of the libraries: ::

             $ LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/lis/lib:/opt/gsl/lib \
                  /path/to/install/location/hophop


