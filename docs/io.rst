Input/Output
============

Read here to learn how to set input parameters to *hophop* and how output is
written.

Input
-----

Parameters are specified in *hophop* on the command line only. The following is
the output of :code:`hophop -h`, that describes all available CLI parameters. ::

    HOP 2.4

    This software simulates hopping in disordered semiconductors with hopping on
    localized states. It uses Monte-Carlo simulation techniques. See the README.rst
    file to learn more.

    Usage: HOP [-h|--help] [-V|--version] [-q|--quiet]
             [-fSTRING|--conf_file=STRING] [-m|--memreq] [--rseed=LONG]
             [-iINT|--nruns=INT] [-P|--parallel] [-tINT|--nthreads=INT]
             [-FFLOAT|--field=FLOAT] [-TFLOAT|--temperature=FLOAT]
             [-lINT|--length=INT] [-XINT|--X=INT] [-YINT|--Y=INT] [-ZINT|--Z=INT]
             [-NINT|--nsites=INT] [-nINT|--ncarriers=INT] [--rc=FLOAT]
             [-pFLOAT|--exponent=FLOAT] [-aFLOAT|--llength=FLOAT] [--gaussian]
             [--lattice] [--removesoftpairs] [--softpairthreshold=FLOAT]
             [--cutoutenergy=FLOAT] [--cutoutwidth=FLOAT]
             [-ILONG|--simulation=LONG] [-RLONG|--relaxation=LONG]
             [-xINT|--nreruns=INT] [--many] [--be] [--mgmres] [--be_it=LONG]
             [--be_oit=LONG] [--tol_abs=FLOAT] [--tol_rel=FLOAT] [--an]
             [-BFLOAT|--percolation_threshold=FLOAT]
             [-oSTRING|--outputfolder=STRING] [--transitions]
             [-ySTRING|--summary=STRING] [-cSTRING|--comment=STRING]

      -h, --help                    Print help and exit
      -V, --version                 Print version and exit
      -q, --quiet                   Don't say anything.  (default=off)
      -f, --conf_file=STRING        Location of a configuration file for the
                                      simulation.
      -m, --memreq                  Estimates the used memory for the specified
                                      parameter set. Print's the information and
                                      exits immediately  (default=off)
          --rseed=LONG              Set the random seed manually.
      -i, --nruns=INT               The number of runs to average over.
                                      (default=`1')
      -P, --parallel                If the runs given with the --nruns option
                                      should be executed using mutliple cores and
                                      parallelization. This suppresses any progress
                                      output of the runs but will be very fast on
                                      multicore systems.  (default=off)
      -t, --nthreads=INT            The number of threads to use during parallel
                                      computing. 0 means all there are.
                                      (default=`0')

    External physical parameters:
      Some external physical quantities.
      -F, --field=FLOAT             The electric field strength in z-direction.
                                      (default=`0.01')
      -T, --temperature=FLOAT       The temperature of the simulation.
                                      (default=`0.3')

    System information:
      Parameters describing the distribution of sites in the system
      -l, --length=INT              This parameter specifies the length of the
                                      (cubic) sample. If it parameter is set, the
                                      options X,Y,Z are ignored!
      -X, --X=INT                   The x-length of the sample. Right now, only
                                      cubic samples should be used, so rather use
                                      the parameter --length.  (default=`50')
      -Y, --Y=INT                   The y-length of the sample. Right now, only
                                      cubic samples should be used, so rather use
                                      the parameter --length.  (default=`50')
      -Z, --Z=INT                   The z-length of the sample. Right now, only
                                      cubic samples should be used, so rather use
                                      the parameter --length.  (default=`50')
      -N, --nsites=INT              The number of localized states. This value has
                                      to be bigger than --ncarriers. Deprecated!
                                      Scale the number os states using --length.
                                      (default=`125000')
      -n, --ncarriers=INT           The number of charge carriers in the system.
                                      (default=`1')
          --rc=FLOAT                Determines up to which distance sites should be
                                      neighbors.  (default=`3')
      -p, --exponent=FLOAT          The exponent of the DOS g(x) = exp(-(x)^p)
                                      (default=`2.0')
      -a, --llength=FLOAT           Localization length of the sites, assumed equal
                                      for all of them.  (default=`0.215')
          --gaussian                Use a Gaussian DOS with std. dev. 1. g(x) =
                                      exp(-1/2*(x)^2)  (default=off)
          --lattice                 Distribute sites on a lattice with distance
                                      unity. Control nearest neighbor hopping and
                                      so on with --rc  (default=off)
          --removesoftpairs         Remove softpairs.  (default=off)
          --softpairthreshold=FLOAT The min hopping rate ratio to define a softpair
                                      (default=`0.95')
          --cutoutenergy=FLOAT      States below this energy will be cut out of the
                                      DOS  (default=`0')
          --cutoutwidth=FLOAT       The width of energies who are cutted.
                                      (default=`0.5')

    Monte carlo simulation:
      The following options matter only, when the system is simulated using a Monte
      Carlo simulation (which is the default)
      -I, --simulation=LONG         The number of hops during which statistics are
                                      collected.  (default=`1000000000')
      -R, --relaxation=LONG         The number of hops to relax.
                                      (default=`100000000')
      -x, --nreruns=INT             How many times should the electron be placed at
                                      some random starting position?  (default=`1')
          --many                    Instead of using the mean field approach,
                                      simulate multiple charge carriers. (slow!!!)
                                      (default=off)

    Balance equations:
      These options only matter, when the solution is found by solving the balance
      equations. (setting the --be flag)
          --be                      Solve balance equations  (default=off)
          --mgmres                  Force use of mgmres instead of lis
                                      (default=off)
          --be_it=LONG              Max inner iterations after which the
                                      calculation is stopped.   (default=`300')
          --be_oit=LONG             Max outer iterations or restarts of the
                                      algorithm.  (default=`10')
          --tol_abs=FLOAT           absolute tolerance for finding the solution
                                      (default=`1e-8')
          --tol_rel=FLOAT           relative tolerance for finding the solution
                                      (default=`1e-8')

    Analytic calculations:
      These options control the analytic calculation of several properties of the
      system, like the transport energy or the mobility.
          --an                      Also try to calculate stuff analytically
                                      (default=off)
      -B, --percolation_threshold=FLOAT
                                    The percolation threshold.  (default=`2.7')

    Output:
      -o, --outputfolder=STRING     The name of the output folder if one wants
                                      output files.
          --transitions             Save all transitions to a file. (Can be big,
                                      scales with -l^3!) Only valid when
                                      --outputfolder is given  (default=off)
      -y, --summary=STRING          The name of the summary file to which one
                                      summary result line is then written.
      -c, --comment=STRING          Specify a string that is appended to the line
                                      in the summary file for better overview over
                                      the simulated data.


Output
------

There are three ways to get output from the simulation:

:code:`-o, --outputfolder`
~~~~~~~~~~~~~~~~~~~~~~~~~~

When the CLI parameter :code:`--outputfolder` (or, equivalently, :code:`-o`) is
specified, *hophop* creates a directory with that value and writes results.

The following files are written:

* :code:`params.conf`:
    A file with the command line parameters given for that simulation. A simulation
    can be started from such a file using the CLI parameter :code:`-f, --conf_file`.
* :code:`1/results.dat`:
    A column-based text file with some simulation parameters and results. Each simulation
    is one line. When the file already exists, a new line will be added. The descriptions
    of the columns are given in the first two lines of the file.

    When multiple runs are simulated, with the parameter :code:`-i, --nruns`, then
    a folder is created for each run, e.g., :code:`1/results.dat`, :code:`2/results.dat`
    etc.
* :code:`1/sites.dat`:
    The generated system and the number of times each site was visited. The columns of the
    file are as follows: ::

        x   y   z   energy    times_visited     times_visited_upward

    :code:`times_visited_upward` is the number of times this site was visited in a hop, where
    the original energy is lower than that of the target site (i.e., the hop is energetically
    an `upward` hop).

    :code:`times_visited` and :code:`times_visited_upward` are only non-zero in KMC mode.

    When multiple runs are simulated, with the parameter :code:`-i, --nruns`, then
    a folder is created for each run, e.g., :code:`1/sites.dat`, :code:`2/sites.dat`
    etc.

:code:`-y, --summary`
~~~~~~~~~~~~~~~~~~~~~

Path to a single columnar summary file, in which system parameters and results are written.
Each simulation is one line. When the file already exists, a new line will be added.
The descriptions of the columns are given in the first two lines of the file.

In BE mode, some columns will be NaN or zero.

stdout
~~~~~~

Some results will also be written to stdout.