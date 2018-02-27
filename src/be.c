/*
 * hophop: Charge transport simulations in disordered systems
 *
 * Copyright (c) 2012-2018 Jan Oliver Oelerich <jan.oliver.oelerich@physik.uni-marburg.de>
 * Copyright (c) 2012-2018 Disordered Many-Particle Physics Group, Philipps-Universität Marburg, Germany
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include "hop.h"
#include "../3rdparty/mgmres/mgmres.h"

#ifdef WITH_LIS

#include "lis.h"
int solve_lis(Site * sites, RunParams * runprms, int nnz, double * x);

#endif /* WITH_LIS */

void BE_solve (Site * sites, Results * res, RunParams * runprms);
int solve_mgmres(Site * sites, RunParams * runprms, int nnz, double * x);

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

    int i, nnz, j, it;
    
    // find out the number of entries
    nnz = 2 * runprms->nSites - 1;
    for (i = 1; i < runprms->nSites; ++i)
        nnz += sites[i].nNeighbors;

    double *x = malloc (sizeof (double) * runprms->nSites);

    if(prms.mgmres)
    {
        it = solve_mgmres(sites, runprms, nnz, x);
    }

#ifdef WITH_LIS

    else
    {
        it = solve_lis(sites, runprms, nnz, x);
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
    free (x);

}

#ifdef WITH_LIS

int
solve_lis(Site * sites, RunParams * runprms, int nnz, double * x)
{
    int i, *ia, *ja, k, it;
    double *a, *rhs;
    SLE *neighbor;

    //allocate memory
    ia = malloc (sizeof (int) * nnz);
    ja = malloc (sizeof (int) * nnz);
    a = malloc (sizeof (double) * nnz);
    rhs = calloc (runprms->nSites, sizeof (double));
    
    // create the arrays
    int counter = 0;

    for (i = 0; i < runprms->nSites; ++i)
    {
        a[counter] = 1;
        ja[counter] = 0;
        ia[counter] = i;
        counter++;

        // diagonal element is always the rate sum
        if(i > 0)
        {
            a[counter] = -1 * sites[i].rateSum;
            ja[counter] = i;
            ia[counter] = i;
            counter++;
        }

        // now the neighbors
        for (k = 0; k < sites[i].nNeighbors; ++k)
        {
            neighbor = &(sites[i].neighbors[k]);

            if(neighbor->s->index == 0)
                continue;
            
            a[counter] = neighbor->rate;
            ja[counter] = neighbor->s->index;
            ia[counter] = i;
            counter++;
        }
    }

    // create the right hand side
    rhs[0] = 1;

    // create the initial guess
    for (i = 0; i < runprms->nSites; ++i)
        x[i] = 1. / runprms->nSites;

    int threads = omp_get_num_threads();
    omp_set_num_threads(1);

    char options[200];
    sprintf(options, 
        "-i gmres -p ilu -tol %e -t %d -restart %d -maxiter %d", 
        prms.be_abs_tol, 1, prms.be_it, prms.be_outer_it * prms.be_it);

    LIS_INT argc = 1;
    char ** argv = (char *[]){""};

    lis_initialize(&argc, &argv);

    LIS_MATRIX lis_A;
    lis_matrix_create(0,&lis_A);
    lis_matrix_set_type(lis_A, LIS_MATRIX_COO);
    lis_matrix_set_size(lis_A,0,runprms->nSites);

    lis_matrix_set_coo(nnz,ja,ia,a,lis_A);
    lis_matrix_assemble(lis_A);

    LIS_VECTOR lis_b, lis_x;
    LIS_SOLVER solver;

    lis_vector_duplicate(lis_A, &lis_b);
    lis_vector_duplicate(lis_A, &lis_x);

    for(i = 0; i < runprms->nSites; ++i)
        lis_vector_set_value(LIS_INS_VALUE,i,rhs[i],lis_b);

    lis_solver_create(&solver);
    lis_solver_set_option(options,solver);
    lis_solve(lis_A,lis_b,lis_x,solver);
    lis_solver_get_iter(solver,&it);

    for(i = 0; i < runprms->nSites; ++i)
        lis_vector_get_value(lis_x, i, &(x[i]));

    lis_finalize();

    omp_set_num_threads(threads);


    free (a);
    free (ia);
    free (ja);
    free (rhs);

    return it;

}

#endif /* WITH_LIS */

int
solve_mgmres(Site * sites, RunParams * runprms, int nnz, double * x)
{
    int i, *ia, *ja, j, k, it;
    double *a, *rhs;
    SLE *neighbor, *neighbor2;

    //allocate memory
    ia = malloc (sizeof (int) * (runprms->nSites + 1));
    ja = malloc (sizeof (int) * nnz);
    a = malloc (sizeof (double) * nnz);
    rhs = calloc (runprms->nSites, sizeof (double));
    
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
    it = pmgmres_ilu_cr (runprms->nSites, nnz, ia, ja, a, x, rhs, 
            prms.be_outer_it, prms.be_it, prms.be_abs_tol, prms.be_rel_tol);


    free (a);
    free (ia);
    free (ja);
    free (rhs);

    return it;
}