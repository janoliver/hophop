#include "hop.h"
#include "mgmres.h"

#include "lis.h"

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
    output (O_SERIAL, "\tSolving balance equations ...");
    fflush (stdout);

    // measure the cpu time
    struct timeval start, end, result;
    gettimeofday (&start, NULL);

    int i, *ia, *ja, nnz, j, k, it;
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

    if(prms.mgmres)
    {
        // perform the calculation
        it = pmgmres_ilu_cr (runprms->nSites, nnz, ia, ja, a, x, rhs, 
            prms.be_outer_it, prms.be_it, prms.be_abs_tol, prms.be_rel_tol);

    }

#ifdef WITH_LIS

    if(!prms.mgmres)
    {
        char options[200];
        sprintf(options, "-i gmres -p ilu -tol %e -t %d -restart %d -maxiter %d", 
            prms.be_abs_tol, 1, prms.be_outer_it, prms.be_it);

        LIS_INT argc = 1;
        char ** argv = (char *[]){""};

        lis_initialize(&argc, &argv);

        LIS_MATRIX lis_A;
        lis_matrix_create(0,&lis_A);
        lis_matrix_set_size(lis_A,0,runprms->nSites);

        lis_matrix_set_csr(nnz,ia,ja,a,lis_A);
        lis_matrix_assemble(lis_A);

        LIS_VECTOR lis_b, lis_x;
        LIS_SOLVER solver;

        lis_vector_duplicate(lis_A, &lis_b);
        lis_vector_duplicate(lis_A, &lis_x);

        for(i = 0; i < runprms->nSites; ++i)
        {
            lis_vector_set_value(LIS_INS_VALUE,i,rhs[i],lis_b);
        }

        lis_solver_create(&solver);
        lis_solver_set_option(options,solver);
        lis_solve(lis_A,lis_b,lis_x,solver);
        lis_solver_get_iters(solver,&it);

        for(i = 0; i < runprms->nSites; ++i)
        {
            lis_vector_get_value(lis_x, i, &(x[i]));
        }

        lis_finalize();
    }

#endif /* WITH_LIS */

    //timer
    gettimeofday (&end, NULL);
    timeval_subtract (&result, &start, &end);
    double elapsed = result.tv_sec + (double) result.tv_usec / 1e6;

    // output
    output (O_SERIAL, "\tDone! %f s duration, %d gmres iterations\n", elapsed, it);
    output (O_PARALLEL, "Finished %d. Iteration (total %d): %fs duration, %d gmres iterations\n",
            runprms->iRun, prms.number_runs, elapsed, it);

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
