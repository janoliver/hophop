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
    Results res;
    MC_simulation (sites, carriers, &res, iRun);
    MC_calculateResults (sites, carriers, &res);

    // add results to the total result struct for averaging
    total[*iRun - 1].mobility = res.mobility;
    total[*iRun - 1].diffusivity = res.diffusivity;
    total[*iRun - 1].currentDensity = res.currentDensity;
    total[*iRun - 1].simulationTime = res.simulationTime;
    total[*iRun - 1].equilibrationEnergy = res.equilibrationEnergy;
    total[*iRun - 1].avgenergy = res.avgenergy;
    total[*iRun - 1].einsteinrelation = res.einsteinrelation;
    
    // write output files
    if (strArgGiven (prms.output_folder))
    {
        writeResults (&res, *iRun);
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
