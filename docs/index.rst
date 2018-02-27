hophop
======

This software simulates charge carrier movement through disordered
systems. It is capable of using the variable range mechanism as well as
some lattice structure. The energetic disorder can be specified as an
input parameter.
The software determines some key features of the charge transport such as
mobility, diffusivity, certain important energies and so on.

Requirements
------------
The following libraries and software is required to compile successfully:

* `gsl <https://www.gnu.org/software/gsl/>`_
    The gnu scientific library provides many mathematical algorithms,
    functions, constants and so on. It is heavily used within the program,
    mainly for pseudo random number generation and probability
    distributions.

* `OpenMP <http://www.openmp.org/>`_
    This is a openmp library for parallel computing. Parts of the software
    are parallelized which is why this is needed. This might become
    optional in the near future.

* `Lis <http://www.ssisc.org/lis/>`_ (optional, but recommended!)
    Lis (Library of Iterative Solvers for linear systems) is used as a solver
    for the balance equations method. When it is not found, the `mgmres` solver
    is used, the source code of which is shipped with `hophop`. However, we
    recommend using Lis where possible.


How to build
------------
Run the following commands to compile the program:
::

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ make install

You may use the :code:`CMake` variables :code:`GSL_ROOT` and :code:`LIS_ROOT`
to specify the directories of GSL and Lis install locations.


Usage
-----
Please see the output of
::

    $ ./hophop -h


Units
-----

* Length
    The parameters specifying the length of the sample, are given in units
    of :math:`N^(-1/3)`. So choosing `-l20` would result in
    :math:`20\times 20\times 20` sites. :math:`N^(-1/3)`
    is thus also the length scale in the simulation. Internally, the sample
    is split into cells of with `2 * --rc * --llength` (cutoff radius *
    loc. length).
* Energies/Temperatures
    All energies are measured in units of the disorder parameter \sigma.
    Since:

    .. math::

       \sigma/kT \approx 3.0

    is a realistic value, the temperature T is usually around :math:`0.3`.


Credits
-------

  * We acknowledge the creators of the supplementary libraries that *hophop* depends on.
  * Many of the algorithms were contributed by Dr. Alexey Nenashev.

.. include:: ../CITATION.rst

.. include:: ../CONTRIBUTING.rst