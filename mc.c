//  Created by Jan Oliver Oelerich

#include "mc.h"

void MC_run(Results * total, int * iRun)
{
    Site    * sites    	   = NULL;
    Carrier * carriers     = NULL;    
    
    // setup random number generator
    RSEED = time(NULL) + *iRun;
    curtime = RSEED;
    if(args.rseed_given)
        RSEED = (unsigned long) args.rseed_arg + *iRun - 1;
    gsl_rng_set(r, RSEED);

    // some output
    if(!args.quiet_given)
        printf("\nRunning %d. iteration (total %d):\n", *iRun, args.nruns_arg);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites();
    MC_createHoppingRates(sites);
    if(args.removesoftpairs_given)
        MC_removeSoftPairs(sites);
    carriers = MC_distributeCarriers(sites);

    // simulate
    Results res;
    MC_simulation(sites, carriers, &res);
    MC_calculateResults(sites, carriers, &res);
    
    // add results to the total result struct for averaging
    total[*iRun-1].mobility        = res.mobility;
    total[*iRun-1].diffusivity     = res.diffusivity;
    total[*iRun-1].fermiEnergy     = res.fermiEnergy;
    total[*iRun-1].transportEnergy = res.transportEnergy;
    total[*iRun-1].currentDensity  = res.currentDensity;
    total[*iRun-1].simulationTime  = res.simulationTime;

    // write output files
    if(args.outputfolder_given)
    {
        writeResults(&res, *iRun);
        writeConfig(*iRun);
        writeSites(sites,*iRun);
        writeSitesConfig(sites,*iRun);
        if(args.transitions_given)
            writeTransitions(sites,*iRun);
    }

    // free resources
    //gsl_rng_free(r);
    SLE * neighbor, * tmp;
    int i;
    for(i = 0; i < args.nsites_arg; ++i)
    {
        // free neighbor memory
        neighbor = sites[i].neighbors;
        while(neighbor)
        {
            tmp = neighbor->next;
            free(neighbor);
            neighbor = tmp;
        }
    }
    free(sites);
    free(carriers);
    //gsl_rng_free(r);
    
    return;
}
