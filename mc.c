//  Created by Jan Oliver Oelerich

#include "hop.h"

void
MC_run (Results * total, int *iRun)
{
    Site *sites = NULL;
    Carrier *carriers = NULL;

    // setup random number generator
    prms.rseed_used = time (NULL) + *iRun;
    prms.curtime = prms.rseed_used;
    if (prms.rseed != 0)
        prms.rseed_used = (unsigned long) prms.rseed + *iRun - 1;
    gsl_rng_set (prms.r, prms.rseed_used);

    // some output
    if (!serialOutput () && !prms.quiet)
        printf ("Starting %d. Iteration (total %d): Thread ID %d\n", *iRun,
                prms.number_runs, omp_get_thread_num ());
    if (serialOutput ())
        printf ("\nRunning %d. iteration (total %d):\n", *iRun,
                prms.number_runs);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites ();
    MC_createHoppingRates (sites);
    if (prms.removesoftpairs)
        MC_removeSoftPairs (sites);
    carriers = MC_distributeCarriers (sites);

    // simulate
    MC_simulation (sites, carriers, &(total[*iRun - 1]), iRun);
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
