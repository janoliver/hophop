//  Created by Jan Oliver Oelerich

#include "hop.h"

// the program parameters
Params prms;

void printResults (Results * results);
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
    Results res;
    init_results (&res);

    // set the number of threads
    omp_set_num_threads (prms.nthreads);

#pragma omp parallel if(prms.parallel) shared(res) private(iRun)
    {
#pragma omp for schedule(dynamic)
        for (iRun = 1; iRun <= prms.number_runs; iRun++)
        {
            // setup random number generator
            RunParams runprms;
            runprms.r = gsl_rng_alloc (prms.T);
            runprms.rseed_used = time (NULL) * iRun;
            if (prms.rseed != 0)
                runprms.rseed_used = (unsigned long) prms.rseed + iRun - 1;
            gsl_rng_set (runprms.r, runprms.rseed_used);
            runprms.nHops = 0;
            runprms.nFailedAttempts = 0;
            runprms.stat = false;
            runprms.simulationTime = 0;
            runprms.iRun = iRun;

            // here is where el magico happens
            if (prms.balance_eq)
            {
                BE_run (&res, &runprms);
            }
            else
            {
                MC_run (&res, &runprms);
            }

            // free the RNG
            gsl_rng_free (runprms.r);

        }
    }

    // average results (see helper.c)
    average_errors (&res);

    // summary
    if (strArgGiven (prms.output_summary))
        writeSummary (&res);

    // output results to the command line
    printResults (&res);

    // free results structs
    //int i;
    //for(i = 0; i< sizeof(res)
    free_results(&res);

    return 0;
}

/*
 * This function prints out a useless header for the program.
 */
void
printApplicationHeader ()
{
    output (O_BOTH, "\n\n################### %s v%s ################### \n",
            PKG_NAME, PKG_VERSION);
}

/*
 * This function just prints out the parameters of the simulation.
 */
void
printSettings ()
{
    // Settings output
    output (O_BOTH, "\nSettings:\n");
    output (O_BOTH, "\t3D sample size: \t\t%d x %d x %d\n",
            prms.length_x, prms.length_y, prms.length_z);
    output (O_BOTH, "\tDOS exponent: \t\t\tp = %1.1f\n", prms.exponent);
    output (O_BOTH, "\tLocalization length of sites: \ta = %2.4f\n",
            prms.loclength);
    output (O_BOTH, "\tTemperature: \t\t\tT = %2.4f\n", prms.temperature);
    output (O_BOTH, "\tField strength: \t\tF = %2.4f\n", prms.field);
    output (O_BOTH, "\tCut-off radius: \t\tr = %2.4f\n", prms.cutoff_radius);

    output (O_PARALLEL, "\tParallelization: \t\tRunning on %d cores\n",
            prms.nthreads);
    output (O_SERIAL, "\tParallelization: \t\tOff\n");
    output (O_BOTH, "\tRealizations for Averaging: \ti = %d\n",
            prms.number_runs);
    output (O_BOTH, "\tMode: \t\t\t\t%s\n\n",
            prms.balance_eq ? "Balance equations" : "Monte carlo simulation");

    if (prms.balance_eq)
    {
        output (O_BOTH, "\tMax. nr. of outer iterations: \t%d\n",
                prms.be_outer_it);
        output (O_BOTH, "\tMax. nr. of inner iterations: \t%d\n", prms.be_it);
        output (O_BOTH, "\tRel. convergence tolerance: \t%e\n",
                prms.be_rel_tol);
        output (O_BOTH, "\tAbs. convergence tolerance: \t%e\n",
                prms.be_abs_tol);
    }
    else
    {
        output (O_BOTH, "\tNumber of reruns:\t\tx = %d\n", prms.number_reruns);
        output (O_BOTH, "\tNumber of carriers: \t\tn = %d\n", prms.ncarriers);
        output (O_BOTH, "\tHops of relaxation: \t\tR = %lu\n", prms.relaxation);
        output (O_BOTH, "\tHops of simulation: \t\tI = %lu\n", prms.simulation);

    }

    output (O_PARALLEL, "\n");
}

/*
 * This function prints out the results of the simulation.
 */
void
printResults (Results * results)
{
    // Results output
    output (O_BOTH, "\nResults:\n");
    output (O_BOTH, "\tMobility in field-direction: \tu   = %e (+- %e)\n",
            results->mobility.avg, results->mobility.err);

    if (prms.balance_eq)
    {
        output (O_BOTH, "\n");
        return;
    }

    output (O_BOTH, "\tDiffusivity perp. to field:  \tD   = %e (+- %e)\n",
            results->diffusivity.avg, results->diffusivity.err);
    output (O_BOTH, "\tEinstein rel. perp. to field: \tu/D = %e (+- %e) e/o\n",
            results->einsteinrelation.avg, results->einsteinrelation.err);
    output (O_BOTH, "\tCurrent density (z-dir):  \tj   = %e (+- %e)\n",
            results->currentDensity.avg, results->currentDensity.err);
    output (O_BOTH, "\tEquilibration Energy: \t\tE_i = %e\n",
            results->equilibrationEnergy.avg);
    output (O_BOTH, "\tSimulated time: \t\tt   = %e\n\n",
            results->simulationTime.avg);

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

    output (O_BOTH, "Estimated memory usage: %5.2f MB\n", mem / (1024 * 1024));
    return;
}
