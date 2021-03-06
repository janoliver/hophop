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

void
MC_run (Results * res, RunParams * runprms)
{
    Site *sites = NULL;
    Carrier *carriers = NULL;
    int i;

    struct timeval start, end, result;

    // some output
    output (O_PARALLEL, "Starting %d. Iteration (total %d): Thread ID %d\n",
            runprms->iRun, prms.number_runs, omp_get_thread_num ());
    output (O_SERIAL, "\nRunning %d. iteration (total %d):\n", runprms->iRun,
            prms.number_runs);

    // create the sites, cells, carriers, hopping rates
    sites = MC_createSites (runprms);
    MC_createHoppingRates (sites, runprms);
    if (prms.removesoftpairs)
        MC_removeSoftPairs (sites, runprms);
    carriers = MC_createCarriers ();

    gettimeofday (&start, NULL);

    // simulate
    for (i = 0; i < prms.number_reruns; ++i)
    {
        MC_distributeCarriers (carriers, sites, runprms);
        MC_simulation (sites, carriers, runprms, i + 1);
    }

    // some more output
    gettimeofday (&end, NULL);
    timeval_subtract (&result, &start, &end);

    double elapsed = result.tv_sec + (double) result.tv_usec / 1e6;

    output (O_PARALLEL,
            "Finished %d. Iteration (total %d): %lu successful hops/sec (%ld failed)\n",
            runprms->iRun, prms.number_runs,
            (size_t) (prms.number_reruns * (prms.relaxation + prms.simulation) /
                      elapsed), runprms->nFailedAttempts);
    output (O_SERIAL, " Done. %lu successful hops/sec (%ld failed)\n",
            (size_t) (prms.number_reruns * (prms.relaxation + prms.simulation) /
                      elapsed), runprms->nFailedAttempts);

    // calculate the results
    MC_calculateResults (sites, carriers, res, runprms);

    // write output files
    if (strArgGiven (prms.output_folder))
    {
        writeResults (res, runprms);
        writeConfig (runprms);
        writeSites (sites, runprms);
        if (prms.output_transitions)
            writeTransitions (sites, runprms);
    }

    // free resources
    for (i = 0; i < runprms->nSites; ++i)
    {
        // free neighbor memory
        free (sites[i].neighbors);
    }
    free (sites);
    free (carriers);

    return;
}
