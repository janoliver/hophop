#include "hop.h"
#include "mgmres.h"

void
BE_run (Results * total, int *iRun)
{
    Site *sites = NULL;

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
        printf ("\nRunning %d. iteration (total %d)\n", *iRun,
                prms.number_runs);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites ();
    MC_createHoppingRates (sites);
    if (prms.removesoftpairs)
        MC_removeSoftPairs (sites);

    // solve
    Results res;
    BE_solve (sites, &res, iRun);

    // add results to the total result struct for averaging
    total[*iRun - 1].mobility = res.mobility;
    total[*iRun - 1].diffusivity = res.diffusivity;
    total[*iRun - 1].currentDensity = res.currentDensity;
    total[*iRun - 1].simulationTime = res.simulationTime;
    total[*iRun - 1].equilibrationEnergy = res.equilibrationEnergy;

    // write output files
    if (strArgGiven (prms.output_folder))
    {
        writeResults (&res, *iRun);
        writeConfig (*iRun);
        writeSites (sites, *iRun);
        writeSitesConfig (sites, *iRun);
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

    return;
}

// preconditioner
void
BE_solve (Site * sites, Results * res, int *iRun)
{
    if (serialOutput ())
    {
        printf ("\tSolving balance equations...");
        fflush (stdout);
    }

    // measure the cpu time
    struct timeval start, end, result;
    struct timezone tz;
    gettimeofday (&start, &tz);

    // create the matrix in sparse triplet form
    int i, *ia, *ja, nnz;
    double *a, *rhs, *x;
    SLE *neighbor, *neighbor2;

    // find out the number of entries
    nnz = 2 * prms.nsites - 1;
    for (i = 1; i < prms.nsites; ++i)
    {
        neighbor = sites[i].neighbors;
        while (neighbor)
        {
            nnz++;
            neighbor = neighbor->next;
        }
    }

    //allocate memory
    ia = malloc (sizeof (int) * (prms.nsites+1));
    ja = malloc (sizeof (int) * nnz);
    a = malloc (sizeof (double) * nnz);
    rhs = calloc (prms.nsites, sizeof (double));
    x = malloc (sizeof (double) * prms.nsites);

    // create the arrays
    int counter = 0;
    ia[0] = 0;
    for(i = 0; i < prms.nsites; ++i)
    {
        // first element is always 1
        a[counter] = 1;
        ja[counter] = i;
        counter++;
    }
    
    for (i = 1; i < prms.nsites; ++i)
    {
        // where does the current row begin?
        ia[i] = counter;

        // diagonal element is always the rate sum
        a[counter] = -1 * sites[i].rateSum;
        ja[counter] = i;
        counter++;
        
        // now the neighbors
        neighbor = sites[i].neighbors;
        while (neighbor)
        {
            // this is an ugly, ugly hack because we need compressed
            // row storage, which is exactly the wrong order compared
            // to how our matrix has to be built. 
            neighbor2 = neighbor->s->neighbors;
            while (neighbor2)
            {
                if(neighbor2->s->index == i) {
                    a[counter] = neighbor2->rate;
                    break;
                }
                neighbor2 = neighbor2->next;
            }
            ja[counter] = neighbor->s->index;
            counter++;
            neighbor = neighbor->next;
        }
    }
    ia[prms.nsites] = counter;

    // create the right hand side
    rhs[0] = 1;

    // create the initial guess
    for (i = 0; i < prms.nsites; ++i)
    {
        x[i] = 1. / prms.nsites;
    }

    // perform the calculation
    pmgmres_ilu_cr (prms.nsites, nnz, ia, ja, a, x, rhs, 1000,
                    1000, prms.be_abs_tol, prms.be_rel_tol);

    //timer
    gettimeofday (&end, &tz);
    timeval_subtract (&result, &start, &end);
    double elapsed = result.tv_sec + (double) result.tv_usec / 1e6;

    // output
    if (serialOutput ())
        printf ("\tDone! %f s duration\n", elapsed);
    if (!serialOutput () && !prms.quiet)
        printf ("Finished %d. Iteration (total %d): %f s duration\n",
                *iRun, prms.number_runs, elapsed);

    // calculate mobility
    double sum = 0;
    for (i = 0; i < prms.nsites; ++i)
    {
        neighbor = sites[i].neighbors;
        while (neighbor)
        {
            sum += x[i] * neighbor->rate * (neighbor->dist.z);
            neighbor = neighbor->next;
        }
    }
    res->mobility = sum / prms.field;

    // free
    free (a);
    free (ia);
    free (ja);
    free (rhs);
    free (x);

}