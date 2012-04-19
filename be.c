#include "hop.h"
#include "mgmres.h"

void
BE_run (Results * total, RunParams * runprms, int *iRun)
{
    Site *sites = NULL;

    // some output
    output (O_PARALLEL, "Starting %d. Iteration (total %d): Thread ID %d\n",
            *iRun, prms.number_runs, omp_get_thread_num ());
    output (O_SERIAL, "\nRunning %d. iteration (total %d)\n", *iRun,
            prms.number_runs);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites (runprms);
    MC_createHoppingRates (sites);
    if (prms.removesoftpairs)
        MC_removeSoftPairs (sites);

    // solve
    BE_solve (sites, &(total[*iRun - 1]), iRun);

    // write output files
    if (strArgGiven (prms.output_folder))
    {
        writeResults (&(total[*iRun - 1]), *iRun);
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
        free (sites[i].neighbors);
    }
    free (sites);

    return;
}

// preconditioner
void
BE_solve (Site * sites, Results * res, int *iRun)
{
    output (O_SERIAL, "\tSolving balance equations...");
    fflush (stdout);

    // measure the cpu time
    struct timeval start, end, result;
    struct timezone tz;
    gettimeofday (&start, &tz);

    // create the matrix in sparse triplet form
    int i, *ia, *ja, nnz, j, k;
    double *a, *rhs, *x;
    SLE *neighbor, *neighbor2;

    // find out the number of entries
    nnz = 2 * prms.nsites - 1;
    for (i = 1; i < prms.nsites; ++i)
        nnz += sites[i].nNeighbors;

    //allocate memory
    ia = malloc (sizeof (int) * (prms.nsites + 1));
    ja = malloc (sizeof (int) * nnz);
    a = malloc (sizeof (double) * nnz);
    rhs = calloc (prms.nsites, sizeof (double));
    x = malloc (sizeof (double) * prms.nsites);

    // create the arrays
    int counter = 0;
    ia[0] = 0;
    for (i = 0; i < prms.nsites; ++i)
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
        for (k = 0; i < sites[k].nNeighbors; ++k)
        {
            neighbor = &(sites[k].neighbors[k]);

            // this is an ugly, ugly hack because we need compressed
            // row storage, which is exactly the wrong order compared
            // to how our matrix has to be built. 
            for (j = 0; j < neighbor->s->nNeighbors; ++j)
            {
                neighbor2 = &(neighbor->s->neighbors[j]);

                if (neighbor2->s->index == i)
                {
                    a[counter] = neighbor2->rate;
                    break;
                }
            }
            ja[counter] = neighbor->s->index;
            counter++;
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
    pmgmres_ilu_cr (prms.nsites, nnz, ia, ja, a, x, rhs, prms.be_outer_it,
                    prms.be_it, prms.be_abs_tol, prms.be_rel_tol);

    //timer
    gettimeofday (&end, &tz);
    timeval_subtract (&result, &start, &end);
    double elapsed = result.tv_sec + (double) result.tv_usec / 1e6;

    // output
    output (O_SERIAL, "\tDone! %f s duration\n", elapsed);
    output (O_PARALLEL, "Finished %d. Iteration (total %d): %f s duration\n",
            *iRun, prms.number_runs, elapsed);

    // calculate mobility
    double sum = 0;
    for (i = 0; i < prms.nsites; ++i)
        for (j = 0; j < sites[i].nNeighbors; ++j)
            sum +=
                x[i] * sites[i].neighbors[j].rate *
                (sites[i].neighbors[j].dist.z);

    res->mobility = sum / prms.field;

    // free
    free (a);
    free (ia);
    free (ja);
    free (rhs);
    free (x);

}
