//  Created by Jan Oliver Oelerich

#include "mc.h"

void hoppingStep(Site * sites, Carrier * carriers, Results * res,
                 bool stat, ALT * tables);
void hop(Carrier * c, SLE * dest, Vector * dist,
         Results * res, Carrier * carriers, bool stat);
void updateCarrier(Carrier * carriers);
int  getTransitionLookupTables(Site * sites, ALT * tables);

/**
   This function runs the iteration of the simulation.  It keeps track of
   the passed "time" (in arbitrary units) and ouputs the progress.
*/
void
MC_simulation(Site * sites, Carrier * carriers, Results * res, int * iRun)
{
    size_t j;
    struct timeval start, end, result;
    struct timezone tz;

    // initialize some results
    res->simulationTime = 0.0;
    res->nHops          = 0;
    res->nFailedAttempts= 0;

    // if the accept/reject method is chosen, initialize tables
    ALT alias_tables;
    if(args.ar_given)
    {
        int n = getTransitionLookupTables(sites, &alias_tables);
        alias_tables.tab = gsl_ran_discrete_preproc(n, alias_tables.weights);
    }

    // relaxation, no time or hop counting
    for(j = 0; j < 100; j++)
    {
        while(res->simulationTime < args.relaxationtime_arg / 100 * j)
            hoppingStep(sites, carriers, res, false, &alias_tables);

        if(!args.quiet_given && (!args.parallel_given || args.nruns_arg == 1))
        {
            printf("\r\tRelaxing...:  \t\t\t%2d%%", (int) j);
            fflush(stdout);
        }
    }
    if(!args.quiet_given && (!args.parallel_given || args.nruns_arg == 1))
        printf(" Done.\n");

    // actual simulation, time and hop counting
    for(j = 0; j < args.ncarriers_arg; ++j)
        carriers[j].occTime -= res->simulationTime;

    gettimeofday(&start, &tz);
    res->simulationTime = 0.0;
    for(j = 0; j < 100; j++)
    {
        while(res->simulationTime < args.simulationtime_arg / 100 * j)
            hoppingStep(sites, carriers, res, true, &alias_tables);

        if(!args.quiet_given && (!args.parallel_given || args.nruns_arg == 1))
        {
            printf("\r\tSimulating...: \t\t\t%2d%%", (int) j);
            fflush(stdout);
        }
    }
    gettimeofday(&end, &tz);
    timeval_subtract(&result, &start, &end);

    // finish occupation time calculations
    for(j = 0; j < args.nsites_arg; ++j)
        if(sites[j].tempOccTime > 0)
            sites[j].totalOccTime
                += res->simulationTime - sites[j].tempOccTime;

    double elapsed = result.tv_sec + (double)result.tv_usec / 1e6;
    if(!args.quiet_given && args.parallel_given && args.nruns_arg > 1)
        printf("Finished %d. Iteration (total %d): \n\t%d successful hops/s (%ld successful, %ld failed)\n",
               *iRun, args.nruns_arg, (int) (res->nHops / elapsed), res->nHops, res->nFailedAttempts);
    if(!args.quiet_given && (!args.parallel_given || args.nruns_arg == 1))
        printf(" Done. %d successful hops/s (%ld successful, %ld failed)\n",
               (int) (res->nHops / elapsed), res->nHops, res->nFailedAttempts);
    
    // free the memory for the alias method
    if(args.ar_given)
    {
        free(alias_tables.dest);
        free(alias_tables.orig);
        free(alias_tables.weights);
        gsl_ran_discrete_free(alias_tables.tab);
    }

}

/**
    This function executes one step of the simulation. It finds the carrier
    or occupied site with the lowest occupation time and defines this as
    the next electron to hop. It then reduces all the other times of all
    other electrons and finds the destination for the hop. Then, the hop()
    function is called with these parameters.
 */
void
hoppingStep(Site * sites, Carrier * carriers,
            Results * res, bool stat, ALT * tables)
{
    Carrier * c = NULL;
    SLE * dest = NULL;

    if(args.ar_given)
    {
        size_t nextTransitionIndex = gsl_ran_discrete(r, tables->tab);
        //size_t nextTransitionIndex = 10;
        res->simulationTime += tables->total;

        c = tables->orig[nextTransitionIndex]->carrier;
        dest = tables->dest[nextTransitionIndex];
        
        // check if jump is accepted
        if(c != NULL &&
           dest->s->carrier == NULL)
        {
            // do the hopping and write some statistics
            if(stat)
                res->nHops++;
            
            hop(c, dest, &dest->dist, res, carriers, stat);
        }
        else if(dest->s->carrier != NULL && c != NULL && stat)
        {
            res->nFailedAttempts++;
            c->nFailedAttempts++;  
        }
        else if(stat)
            res->nFailedAttempts++;
           
    }
    else
    {
        double randomHopProb, probSum;
    

        // heap
        c = &carriers[0];

        // determine the next destination site
        randomHopProb  = (float) gsl_rng_uniform(r) * c->site->rateSum;
        probSum        = 0.0;
        SLE * neighbor = c->site->neighbors;
        while(neighbor)
        {
            if(probSum > randomHopProb)
                break;

            dest     = neighbor;
            probSum += neighbor->rate;
            neighbor = neighbor->next;
        }

        res->simulationTime = c->occTime;

        if(dest->s->carrier == NULL)
        {

            // do the hopping and write some statistics
            if(stat)
                res->nHops++;

            hop(c, dest, &dest->dist, res, carriers, stat);
        }
        else
        {
            if(stat)
            {
                res->nFailedAttempts++;
                c->nFailedAttempts++;
            }
        }
        updateCarrier(carriers);
    }
}

/**
    In this function, one hop is executed. Here, the track-keeping of
    energies and so on should happen. It updates both sites that take part
    in the hopping process and all of the sites around these two sites.
 */
void
hop(Carrier * c, SLE * dest, Vector * dist,
    Results * res, Carrier * carriers, bool stat)
{
    Site * orig = c->site;

    // update origin site
    orig->carrier = NULL;

    if(stat)
    {
        dest->nTransitions++;
        
        orig->totalOccTime += res->simulationTime - orig->tempOccTime;
        orig->tempOccTime   = 0.0;
        dest->s->tempOccTime   = res->simulationTime;

        c->dx += dist->x;
        c->dy += dist->y;
        c->dz += dist->z;
    }

    // update carrier
    c->site    = dest->s;

    // update destination site
    dest->s->carrier = c;
    if(stat && orig->energy < dest->s->energy)
        dest->s->visitedUpward++;
    else if(stat)
        dest->s->visited++;

}

/**
   Subtract the `struct timeval' values X and Y, storing the result in
   RESULT.  Return 1 if the difference is negative, otherwise 0.
*/
int
timeval_subtract(struct timeval * result,
                 struct timeval * y, struct timeval * x)
{
    /* Perform the carry for the later subtraction by updating y. */
    if (x->tv_usec < y->tv_usec) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000) {
        int nsec = (x->tv_usec - y->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    /* Compute the time remaining to wait.
       tv_usec is certainly positive. */
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;

    /* Return 1 if result is negative. */
    return x->tv_sec < y->tv_sec;
}

/**
   This is a binary tree search for the carrier with the lowest occupation
   time. It moves it to the top of the carriers list. The last jumped
   carrier is assigned a new occupation time.
 */
void
updateCarrier(Carrier * c)
{
    int smallest, i = 0;

    c[0].occTime += (float) gsl_ran_exponential(r, 1.0) / c[0].site->rateSum;

    // children always sit at c[2i+1] and c[2i+2]
    Carrier tmp;
    while(1)
    {
    	smallest = (2*i+1 < args.ncarriers_arg &&
                    c[i].occTime > c[2*i+1].occTime) ? 2*i+1 : i;
    	if(2*i+2 < args.ncarriers_arg &&
    	   c[smallest].occTime > c[2*i+2].occTime)
            smallest = 2*i+2;

    	if(i == smallest)
            break;

    	// swap
    	tmp = c[i];
    	c[i] = c[smallest];
    	c[smallest] = tmp;
    	i = smallest;
    }
}


/**
   Generates the lookup table used for the accept/reject method
   resp. Walker's random number generator.
 */
int
getTransitionLookupTables(Site * sites, ALT * tables)
{
    size_t k = 0, l = 0;
    tables->weights = NULL;
    tables->orig = NULL;
    tables->dest = NULL;
    SLE * neighbor;
    tables->total = 0;
    
    // find memory to be alloced
    for(k = 0; k < args.nsites_arg; ++k)
        l += sites[k].nNeighbors;
   
    tables->weights = malloc(sizeof(double) * l);
    tables->orig    = malloc(sizeof(Site *) * l);
    tables->dest    = malloc(sizeof(SLE *) * l);

    l = 0;
    for(k = 0; k < args.nsites_arg; ++k)
    {
        printf("\r\tBuilding Tables...:  \t\t%2d%%", (int) (100*k/args.nsites_arg));
        fflush(stdout);
                
        if(sites[k].nNeighbors <= 0)
            continue;
        
        neighbor = sites[k].neighbors;
        while(neighbor)
        {
            tables->weights[l] = neighbor->rate;
            tables->orig[l]    = &sites[k];
            tables->dest[l]    = neighbor;

            tables->total += neighbor->rate;
            
            ++l;
            neighbor = neighbor->next;
            
        }
    }
    tables->total = 1./tables->total;
    
    printf(" Done.\n");
    return l;
}
