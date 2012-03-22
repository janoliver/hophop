//  Created by Jan Oliver Oelerich

#include "hop.h"

// the program parameters
Params prms;

void printResults (Results * results, Results * errors);
void printApplicationHeader ();
void printSettings ();
void checkApplicationSettings (int argc, char **argv);
void printEstimatedMemory ();

int
main (int argc, char **argv)
{
    // parse and check command line arguments
    // and generate the params struct which
    // is used all over the software
    generateParams (&prms, argc, argv);

    if (prms.memreq)
    {
        printEstimatedMemory ();
        return 0;
    }

    //output
    printApplicationHeader ();
    printSettings ();

    // start simulation
    int iRun;
    Results total[prms.number_runs], res, error;
    init_results (&res);
    init_results (&error);

    // set the number of threads
    omp_set_num_threads (prms.nthreads);

#pragma omp parallel if(prms.parallel) shared(total) private(iRun)
    {
#pragma omp for schedule(dynamic)
        for (iRun = 1; iRun <= prms.number_runs; iRun++)
        {

            init_results (&(total[iRun - 1]));

            // here is where el magico happens
            if (prms.balance_eq)
            {
                BE_run (total, &iRun);
            }
            else
            {
                MC_run (total, &iRun);
            }

            // always run analytic calculations as well
            AL_run (total, &iRun);
        }
    }

    // average results (see helper.c)
    average_results (&res, total, &error);

    // summary
    if (strArgGiven (prms.output_summary))
        writeSummary (&res, &error);

    printf ("%e", res.analytic_mobility);

    // output results to the command line
    printResults (&res, &error);

    gsl_rng_free (prms.r);

    return 0;
}

/*
 * This function prints out a useless header for the program.
 */
void
printApplicationHeader ()
{
    if (prms.quiet)
        return;
    printf ("\n\n################### %s v%s ################### \n",
            PKG_NAME, PKG_VERSION);
}

/*
 * This function just prints out the parameters of the simulation.
 */
void
printSettings ()
{
    if (prms.quiet)
        return;
    // Settings output
    printf ("\nSettings:\n");
    printf ("\t3D sample size: \t\t%d x %d x %d\n",
            prms.length_x, prms.length_y, prms.length_z);
    printf ("\tDOS exponent: \t\t\tp = %1.1f\n", prms.exponent);
    printf ("\tDOS width: \t\t\ts = %2.4f\n", prms.sigma);
    printf ("\tLocalization length of sites: \ta = %2.4f\n", prms.loclength);
    printf ("\tTemperature: \t\t\tT = %2.4f\n", prms.temperature);
    printf ("\tField strength: \t\tF = %2.4f\n", prms.field);
    printf ("\tCut-off radius: \t\tr = %2.4f\n", prms.cutoff_radius);
    printf ("\tNumber of sites: \t\tN = %d\n\n", prms.nsites);

    if (prms.parallel)
        printf ("\tParallelization: \t\tRunning on %d cores\n", prms.nthreads);
    else
        printf ("\tParallelization: \t\tOff\n");


    printf ("\tMode: \t\t\t\t%s\n",
            prms.balance_eq ? "Balance equations" : "Monte carlo simulation");

    if (prms.balance_eq)
    {
        printf ("\tMax. nr. of outer iterations: \t%d\n", prms.be_outer_it);
        printf ("\tMax. nr. of inner iterations: \t%d\n", prms.be_it);
        printf ("\tRel. convergence tolerance: \t%e\n", prms.be_rel_tol);
        printf ("\tAbs. convergence tolerance: \t%e\n", prms.be_abs_tol);
    }
    else
    {
        printf ("\tNumber of carriers: \t\tn = %d\n", prms.ncarriers);
        printf ("\tHops of relaxation: \t\tR = %lu\n", prms.relaxation);
        printf ("\tHops of simulation: \t\tI = %lu\n", prms.simulation);
        printf ("\tTechnique: \t\t\t%s\n",
                prms.accept_reject ? "Accept/Reject" : "Standard");
    }

    if (!serialOutput () && !prms.quiet)
        printf ("\n");
}

/*
 * This function prints out the results of the simulation.
 */
void
printResults (Results * results, Results * error)
{
    if (prms.quiet)
        return;

    // Results output
    printf ("\nResults:\n");
    printf ("\tMobility in field-direction: \tu   = %e (+- %e)\n",
            results->mobility, error->mobility);

    if (prms.balance_eq)
    {
        printf ("\n");
        return;
    }

    printf ("\tDiffusivity perp. to field:  \tD   = %e (+- %e)\n",
            results->diffusivity, error->diffusivity);
    printf ("\tEinstein rel. perp. to field: \tu/D = %e (+- %e) e/o\n",
            results->einsteinrelation, error->einsteinrelation);
    printf ("\tCurrent density (z-dir):  \tj   = %e (+- %e)\n",
            results->currentDensity, error->currentDensity);
    printf ("\tEquilibration Energy: \t\tE_i = %e\n",
            results->equilibrationEnergy);
    printf ("\tSimulated time: \t\tt   = %e\n\n", results->simulationTime);

}

void
printEstimatedMemory ()
{
    double mem = 0;

    // sites
    mem += prms.nsites * sizeof (Site);

    // carriers
    mem += prms.ncarriers * sizeof (Carrier);

    // neighbor lists
    mem +=
        pow (prms.cutoff_radius,
             3) * 4. / 3. * M_PI * prms.nsites * sizeof (SLE);

    // parallelization
    if (prms.parallel && prms.number_runs >= omp_get_max_threads ())
        mem *= omp_get_max_threads ();
    else if (prms.parallel)
        mem *= prms.number_runs;

    // results
    mem += (prms.number_runs + 2) * sizeof (Results);

    // some offset that shouldn't scale with the system...
    mem += 100 * 1024 * 1024;

    printf ("Estimated memory usage: %5.2f MB\n", mem / (1024 * 1024));
    return;
}
