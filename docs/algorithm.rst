Algorithm
=========

*hophop* simulates hopping charge transport through a 3D system of 
localized states, which we're going to refer to as `sites`. Sites in
*hophop* have no spatial extend.

Charge carriers (electrons or holes) move via incoherent tunnelling transitions
between the sites. Such transitions are called **hop** (hence, the name *hophop*).
The rate for such a transition (between sites :math:`i` and :math:`j` is given
by the Miller-Abrahams expression:

.. math::

    \nu_{ij} = \nu_0 \exp\left\{-\frac{2d_{ij}}{\alpha}\right\}
        \exp\left\{-\frac{\varepsilon_j - \varepsilon_i + |\varepsilon_j - \varepsilon_i|}{2kT}\right\}

where :math:`\nu_0` is of the order of the vibrational frequency of the atoms,
:math:`d_{ij}` is the spatial distance between the sites, :math:`\varepsilon_i` and
:math:`\varepsilon_j` are the sites' energies and :math:`kT` is thermal energy.

In *hophop*, there are two main algorithms that can be used for simulating transport:

* Kinetic Monte Carlo (KMC)
    Highly optimized KMC implementation that is able to simulate multiple charge carriers.
    The memory footprint of the KMC mode is smaller than that of the Balance Equations,
    however, especially for small systems (up to about 1e6 sites), KMC may be perform worse.
    Before statistics about hopping can be collected, one should do a number of relaxation
    transitions so that the system can reach thermal equilibrium.
* Balance Equation approach (BE)
    Solving the linearized BE for the system results in the occupations of each sites in thermal
    equilibrium. The current implementation always assumes an empty system, i.e., a single charge
    carrier. A steady state (thermal equilibrium) is ensured.

For a description of the two algorithms and deeper insights into the theory of
hopping transport, please have a look at **Part I** of
`Jan Oliver Oelerich's PhD Thesis <https://www.staff.uni-marburg.de/~oelericj/theses/Oelerich_PhD.pdf>`_,
in particular **Chapter 7** for a description of the numerics.

Units
-----

*hophop* uses the following units in input and output:

* Length
    The parameters specifying the length of the sample, are given in units
    of :math:`N^{-\frac{1}{3}}`, where :math:`N` is the total number of sites in the
    system. Choosing :code:`-l20` on the command line, for example, would result in
    :math:`20\times 20\times 20` sites. :math:`N^{-\frac{1}{3}}`
    is thus also the length scale in the simulation. Internally, the sample
    is split into cells of with `2 * --rc * --llength` (cutoff radius *
    loc. length).

* Energies/Temperatures
    All energies are measured in units of the disorder parameter :math:`\sigma`.
    Since:

    .. math::

       \sigma/kT \approx 3.0

    is a realistic value, the temperature T is usually around :math:`0.3`.

* Times
    Times are measured in :math:`\nu_0^{-1}`, where :math:`\nu_0` is the prefactor
    of the Miller-Abrahams hopping rates (see the top of this page).

* Charges
    Units of the elementary charge :math:`e`.

All other units are derived from these:

* Electric field:
    :math:`\sigma/eN^{-\frac{1}{3}}`

* Mobility:
    :math:`N^{-\frac{2}{3}}/(e\sigma \nu_0^{-1})`

* and so on...
    For other quantities, please just express them in terms of the above units.