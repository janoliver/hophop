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
average_errors (Results * res)
{
    Result *results[] = {
        &(res->mobility),
        &(res->diffusivity),
        &(res->einsteinrelation),
        &(res->currentDensity),
        &(res->equilibrationEnergy),
        &(res->avgenergy),

        &(res->nHops),
        &(res->nFailedAttempts),
        &(res->simulationTime),
        &(res->nSites)
    };

    int nResults = sizeof (results) / sizeof (Result *);

    int i, j, c;

    // calculate averages
    for (i = 0; i < nResults; ++i)
    {
        c = 0;
        for (j = 0; j < prms.number_runs; ++j)
        {
            if (results[i]->done[j])
            {

                results[i]->avg += results[i]->values[j];

                c++;
            }
        }
        results[i]->avg /= c;
    }

    // calculate errors
    for (i = 0; i < nResults; ++i)
    {
        c = 0;
        for (j = 0; j < prms.number_runs; ++j)
        {
            if (results[i]->done[j])
            {
                results[i]->err +=
                    pow (results[i]->avg - results[i]->values[j], 2);
                c++;
            }
        }
        results[i]->err = sqrt (results[i]->err / c);
    }

    // set finish time
    time(&(res->time_finished));
}

void
init_results (Results * res)
{
    Result *results[] = {
        &(res->mobility),
        &(res->diffusivity),
        &(res->einsteinrelation),
        &(res->currentDensity),
        &(res->equilibrationEnergy),
        &(res->avgenergy),

        &(res->nHops),
        &(res->nFailedAttempts),
        &(res->simulationTime),
        &(res->nSites)
    };

    int nResults = sizeof (results) / sizeof (Result *);

    int i, j;
    for (i = 0; i < nResults; ++i)
    {
        results[i]->avg = 0;
        results[i]->err = 0;
        results[i]->done = (bool *) malloc (sizeof (bool) * prms.number_runs);
        results[i]->values =
            (double *) malloc (sizeof (double) * prms.number_runs);

        for (j = 0; j < prms.number_runs; ++j)
        {
            results[i]->values[j] = 0;
            results[i]->done[j] = false;
        }
    }

    time(&(res->time_start));
}

void
free_results (Results * res)
{

    Result *results[] = {
        &(res->mobility),
        &(res->diffusivity),
        &(res->einsteinrelation),
        &(res->currentDensity),
        &(res->equilibrationEnergy),
        &(res->avgenergy),

        &(res->nHops),
        &(res->nFailedAttempts),
        &(res->simulationTime),
        &(res->nSites)
    };

    int nResults = sizeof (results) / sizeof (Result *);


    int i;
    for (i = 0; i < nResults; ++i)
    {
        free (results[i]->done);
        free (results[i]->values);
    }
}

int
output (int mode, const char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);

    if (prms.quiet && mode != O_FORCE)
        return 0;

    if ((mode == O_FORCE || mode == O_BOTH) ||
        (mode == O_PARALLEL && prms.parallel && prms.number_runs > 1) ||
        (mode == O_SERIAL && (!prms.parallel || prms.number_runs == 1)))
        return vprintf (fmt, args);

    return 0;
}
