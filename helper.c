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
