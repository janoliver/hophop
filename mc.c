//  Created by Jan Oliver Oelerich

#include "mc.h"

void
MC_run (Results * total, int *iRun)
{
    Site *sites = NULL;
    Carrier *carriers = NULL;

    // setup random number generator
    RSEED = time (NULL) + *iRun;
    curtime = RSEED;
    if (args.rseed_given)
        RSEED = (unsigned long) args.rseed_arg + *iRun - 1;
    gsl_rng_set (r, RSEED);

    // some output
    if (!args.quiet_given && args.parallel_given && args.nruns_arg > 1)
        printf ("Starting %d. Iteration (total %d): Thread ID %d\n", *iRun,
                args.nruns_arg, omp_get_thread_num ());
    if (!args.quiet_given && (!args.parallel_given || args.nruns_arg == 1))
        printf ("\nRunning %d. iteration (total %d):\n", *iRun, args.nruns_arg);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites ();
    MC_createHoppingRates (sites);
    if (args.removesoftpairs_given)
        MC_removeSoftPairs (sites);
    carriers = MC_distributeCarriers (sites);

    // simulate
    Results res;
    MC_simulation (sites, carriers, &res, iRun);
    MC_calculateResults (sites, carriers, &res);

    // add results to the total result struct for averaging
    total[*iRun - 1].mobility = res.mobility;
    total[*iRun - 1].diffusivity = res.diffusivity;
    total[*iRun - 1].currentDensity = res.currentDensity;
    total[*iRun - 1].simulationTime = res.simulationTime;
    total[*iRun - 1].equilibrationEnergy = res.equilibrationEnergy;

    // write output files
    if (args.outputfolder_given)
    {
        writeResults (&res, *iRun);
        writeConfig (*iRun);
        writeSites (sites, *iRun);
        writeSitesConfig (sites, *iRun);
        if (args.transitions_given)
            writeTransitions (sites, *iRun);
    }

    // free resources
    SLE *neighbor, *tmp;
    int i;
    for (i = 0; i < args.nsites_arg; ++i)
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
