#include "hop.h"
#include "mgmres.h"

void BE_solve (Site * sites, Results * res, RunParams * runprms);

void
BE_run (Results * res, RunParams * runprms)
{
    Site *sites = NULL;

    // some output
    output (O_PARALLEL, "Starting %d. Iteration (total %d): Thread ID %d\n",
            runprms->iRun, prms.number_runs, omp_get_thread_num ());
    output (O_SERIAL, "\nRunning %d. iteration (total %d)\n", runprms->iRun,
            prms.number_runs);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites (runprms);
    MC_createHoppingRates (sites, runprms);
    if (prms.removesoftpairs)
        MC_removeSoftPairs (sites, runprms);

    // solve
    BE_solve (sites, res, runprms);

    // write output files
    if (strArgGiven (prms.output_folder))
    {
        writeResults (res, runprms);
        writeConfig (runprms);
        writeSites (sites, runprms);
    }

    // free resources
    int i;
    for (i = 0; i < runprms->nSites; ++i)
        free (sites[i].neighbors);
    free (sites);

    return;
}

// preconditioner
void
BE_solve (Site * sites, Results * res, RunParams * runprms)
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
    nnz = 2 * runprms->nSites - 1;
    for (i = 1; i < runprms->nSites; ++i)
        nnz += sites[i].nNeighbors;

    //allocate memory
    ia = malloc (sizeof (int) * (runprms->nSites + 1));
    ja = malloc (sizeof (int) * nnz);
    a = malloc (sizeof (double) * nnz);
    rhs = calloc (runprms->nSites, sizeof (double));
    x = malloc (sizeof (double) * runprms->nSites);

    // create the arrays
    int counter = 0;
    ia[0] = 0;
    for (i = 0; i < runprms->nSites; ++i)
    {
        // first element is always 1
        a[counter] = 1;
        ja[counter] = i;
        counter++;
    }

    for (i = 1; i < runprms->nSites; ++i)
    {
        // where does the current row begin?
        ia[i] = counter;

        // diagonal element is always the rate sum
        a[counter] = -1 * sites[i].rateSum;
        ja[counter] = i;
        counter++;

        // now the neighbors
        for (k = 0; k < sites[i].nNeighbors; ++k)
        {
            neighbor = &(sites[i].neighbors[k]);

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
    ia[runprms->nSites] = counter;

    // create the right hand side
    rhs[0] = 1;

    // create the initial guess
    for (i = 0; i < runprms->nSites; ++i)
    {
        x[i] = 1. / runprms->nSites;
    }

    // perform the calculation
    pmgmres_ilu_cr (runprms->nSites, nnz, ia, ja, a, x, rhs, prms.be_outer_it,
                    prms.be_it, prms.be_abs_tol, prms.be_rel_tol);

    //timer
    gettimeofday (&end, &tz);
    timeval_subtract (&result, &start, &end);
    double elapsed = result.tv_sec + (double) result.tv_usec / 1e6;

    // output
    output (O_SERIAL, "\tDone! %f s duration\n", elapsed);
    output (O_PARALLEL, "Finished %d. Iteration (total %d): %f s duration\n",
            runprms->iRun, prms.number_runs, elapsed);

    // calculate mobility
    double sum = 0;
    for (i = 0; i < runprms->nSites; ++i)
        for (j = 0; j < sites[i].nNeighbors; ++j)
            sum +=
                x[i] * sites[i].neighbors[j].rate *
                (sites[i].neighbors[j].dist.z);

    res->mobility.values[runprms->iRun - 1] = sum / prms.field;
    res->mobility.done[runprms->iRun - 1] = true;

    res->nSites.values[runprms->iRun - 1] = (float)runprms->nSites;
    res->nSites.done[runprms->iRun - 1] = true;


    // free
    free (a);
    free (ia);
    free (ja);
    free (rhs);
    free (x);

}
