//  Created by Jan Oliver Oelerich

#include "hop.h"

void
MC_run (Results * total, RunParams * runprms, int *iRun)
{
    Site *sites = NULL;
    Carrier *carriers = NULL;

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
    carriers = MC_distributeCarriers (sites, runprms);

    // simulate
    MC_simulation (sites, carriers, &(total[*iRun - 1]), runprms, iRun);
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
    int i;
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
