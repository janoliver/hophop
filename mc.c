//  Created by Jan Oliver Oelerich

#include "hop.h"

void
MC_run (Results * total, RunParams * runprms, int *iRun)
{
    Site *sites = NULL;
    Carrier *carriers = NULL;
    int i;

    struct timeval start, end, result;
    struct timezone tz;

    // some output
    output (O_PARALLEL, "Starting %d. Iteration (total %d): Thread ID %d\n",
            *iRun, prms.number_runs, omp_get_thread_num ());
    output (O_SERIAL, "\nRunning %d. iteration (total %d):\n", *iRun,
            prms.number_runs);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites (runprms);
    MC_createHoppingRates (sites);
    if (prms.removesoftpairs)
        MC_removeSoftPairs (sites);
    carriers = MC_createCarriers (sites);

    gettimeofday (&start, &tz);

    // simulate
    for (i = 0; i < prms.number_reruns; ++i)
    {
        MC_distributeCarriers (carriers, sites, runprms, &(total[*iRun - 1]));
        MC_simulation (sites, carriers, &(total[*iRun - 1]), runprms, iRun,
                       i + 1);
    }

    // some more output
    gettimeofday (&end, &tz);
    timeval_subtract (&result, &start, &end);

    double elapsed = result.tv_sec + (double) result.tv_usec / 1e6;

    output (O_PARALLEL,
            "Finished %d. Iteration (total %d): %d successful hops/sec (%ld failed)\n",
            *iRun, prms.number_runs,
            (int) (prms.number_reruns * (prms.relaxation + prms.simulation) /
                   elapsed), total[*iRun - 1].nFailedAttempts);
    output (O_SERIAL, " Done. %d successful hops/sec (%ld failed)\n",
            (int) (prms.number_reruns * (prms.relaxation + prms.simulation) /
                   elapsed), total[*iRun - 1].nFailedAttempts);

    // calculate the results
    MC_calculateResults (sites, carriers, &(total[*iRun - 1]));

    // write output files
    if (strArgGiven (prms.output_folder))
    {
        writeResults (&(total[*iRun - 1]), *iRun);
        writeConfig (*iRun);
        writeSites (sites, *iRun);
        writeSitesConfig (sites, *iRun);
        if (prms.output_transitions)
            writeTransitions (sites, *iRun);
    }

    // free resources
    for (i = 0; i < prms.nsites; ++i)
    {
        // free neighbor memory
        free (sites[i].neighbors);
    }
    free (sites);
    free (carriers);

    return;
}
