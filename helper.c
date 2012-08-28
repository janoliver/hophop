#include "hop.h"

void
average_results (Results * res, Results * total, Results * error)
{
    int iRun;

    // perform the averaging
    for (iRun = 0; iRun < prms.number_runs; iRun++)
    {
        add_results_to_results (res, &(total[iRun]));
    }

    divide_results_by_scalar (res, (double) prms.number_runs);

    // find the errors
    for (iRun = 0; iRun < prms.number_runs; iRun++)
    {
        error->mobility += pow (res->mobility - total[iRun].mobility, 2);
        error->diffusivity +=
            pow (res->diffusivity - total[iRun].diffusivity, 2);
        error->simulationTime +=
            pow (res->simulationTime - total[iRun].simulationTime, 2);
        error->equilibrationEnergy +=
            pow (res->equilibrationEnergy - total[iRun].equilibrationEnergy, 2);
        error->avgenergy += pow (res->avgenergy - total[iRun].avgenergy, 2);
        error->einsteinrelation += pow (res->einsteinrelation -
                                        total[iRun].einsteinrelation, 2);

    }

    error->mobility = sqrt (error->mobility / prms.number_runs);
    error->diffusivity = sqrt (error->diffusivity / prms.number_runs);
    error->currentDensity = sqrt (error->currentDensity / prms.number_runs);
    error->simulationTime = sqrt (error->simulationTime / prms.number_runs);
    error->equilibrationEnergy =
        sqrt (error->equilibrationEnergy / prms.number_runs);
    error->avgenergy = sqrt (error->avgenergy / prms.number_runs);
    error->einsteinrelation = sqrt (error->einsteinrelation / prms.number_runs);

}

void
init_results (Results * res)
{
    res->mobility = 0.0;
    res->diffusivity = 0.0;
    res->currentDensity = 0.0;
    res->simulationTime = 0.0;
    res->equilibrationEnergy = 0.0;
    res->avgenergy = 0.0;
    res->einsteinrelation = 0.0;

    res->nFailedAttempts = 0;
}

void
add_results_to_results (Results * res, Results * to_add)
{
    res->mobility += to_add->mobility;
    res->diffusivity += to_add->diffusivity;
    res->currentDensity += to_add->currentDensity;
    res->simulationTime += to_add->simulationTime;
    res->equilibrationEnergy += to_add->equilibrationEnergy;
    res->avgenergy += to_add->avgenergy;
    res->einsteinrelation += to_add->einsteinrelation;

}

void
divide_results_by_scalar (Results * res, double scalar)
{

    res->mobility /= scalar;
    res->diffusivity /= scalar;
    res->currentDensity /= scalar;
    res->simulationTime /= scalar;
    res->equilibrationEnergy /= scalar;
    res->avgenergy /= scalar;
    res->einsteinrelation /= scalar;

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
