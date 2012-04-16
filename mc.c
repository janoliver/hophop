//  Created by Jan Oliver Oelerich

#include "hop.h"

void
MC_run (Results * total, RunParams * runprms, int *iRun)
{
    Site *sites = NULL;
    Carrier *carriers = NULL;
    int i;
    size_t j;
    struct timeval start, end, result;
    struct timezone tz;

    gettimeofday (&start, &tz);

    // some output
    if (!serialOutput () && !prms.quiet)
        printf ("Starting %d. Iteration (total %d): Thread ID %d\n", *iRun,
                prms.number_runs, omp_get_thread_num ());
    if (serialOutput ())
        printf ("\nRunning %d. iteration (total %d):\n", *iRun,
                prms.number_runs);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites (runprms);
    MC_createHoppingRates (sites);
    if (prms.removesoftpairs)
        MC_removeSoftPairs (sites);
    carriers = MC_createCarriers (sites);
    
    // simulate
    for(i = 0; i < prms.number_reruns; ++i)
    {
        MC_distributeCarriers (carriers, sites, runprms);
        MC_simulation (sites, carriers, &(total[*iRun - 1]), runprms, iRun, i+1);
        total[*iRun - 1].totalSimulationTime += total[*iRun - 1].simulationTime;
    }

    // finish statistics
    for (j = 0; j < prms.nsites; ++j)
        if (sites[j].tempOccTime > 0)
            sites[j].totalOccTime += total[*iRun - 1].totalSimulationTime - sites[j].tempOccTime;
    
    // some more output
    gettimeofday (&end, &tz);
    timeval_subtract (&result, &start, &end);

    double elapsed = result.tv_sec + (double) result.tv_usec / 1e6;
    if (!serialOutput () && !prms.quiet)
        printf
            ("Finished %d. Iteration (total %d): %d successful hops/sec (%ld failed)\n",
             *iRun, prms.number_runs, (int) (prms.number_reruns * (prms.relaxation + prms.simulation) / elapsed),
             total[*iRun - 1].nFailedAttempts);
    if (serialOutput ())
        printf (" Done. %d successful hops/sec (%ld failed)\n",
                (int) (prms.number_reruns * (prms.relaxation + prms.simulation) / elapsed), total[*iRun - 1].nFailedAttempts);

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
    SLE *neighbor, *tmp;
    for (i = 0; i < prms.nsites; ++i)
    {
        // free neighbor memory
        neighbor = sites[i].neighbors;
        while (neighbor)
        {
            tmp = neighbor->next;
            free (neighbor);
            neighbor = tmp;
        }
    }
    free (sites);
    free (carriers);

    return;
}
