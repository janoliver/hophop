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
        error->totalSimulationTime +=
            pow (res->totalSimulationTime - total[iRun].totalSimulationTime, 2);
        error->equilibrationEnergy +=
            pow (res->equilibrationEnergy - total[iRun].equilibrationEnergy, 2);
        error->avgenergy += pow (res->avgenergy - total[iRun].avgenergy, 2);
        error->einsteinrelation += pow (res->einsteinrelation -
                                        total[iRun].einsteinrelation, 2);

        error->analytic_mobility += pow (res->analytic_mobility -
                                         total[iRun].analytic_mobility, 2);
        error->analytic_fermienergy += pow (res->analytic_fermienergy -
                                            total[iRun].analytic_fermienergy,
                                            2);
        error->analytic_transportenergy +=
            pow (res->analytic_transportenergy -
                 total[iRun].analytic_transportenergy, 2);
    }

    error->mobility = sqrt (error->mobility / prms.number_runs);
    error->diffusivity = sqrt (error->diffusivity / prms.number_runs);
    error->currentDensity = sqrt (error->currentDensity / prms.number_runs);
    error->simulationTime = sqrt (error->simulationTime / prms.number_runs);
    error->totalSimulationTime = sqrt (error->totalSimulationTime / prms.number_runs);
    error->equilibrationEnergy =
        sqrt (error->equilibrationEnergy / prms.number_runs);
    error->avgenergy = sqrt (error->avgenergy / prms.number_runs);
    error->einsteinrelation = sqrt (error->einsteinrelation / prms.number_runs);

    error->analytic_mobility =
        sqrt (error->analytic_mobility / prms.number_runs);
    error->analytic_fermienergy =
        sqrt (error->analytic_fermienergy / prms.number_runs);
    error->analytic_transportenergy =
        sqrt (error->analytic_transportenergy / prms.number_runs);

}

void
init_results (Results * res)
{
    res->mobility = 0.0;
    res->diffusivity = 0.0;
    res->currentDensity = 0.0;
    res->simulationTime = 0.0;
    res->totalSimulationTime = 0.0;
    res->equilibrationEnergy = 0.0;
    res->avgenergy = 0.0;
    res->einsteinrelation = 0.0;

    res->analytic_mobility = 0.0;
    res->analytic_fermienergy = 0.0;
    res->analytic_transportenergy = 0.0;
}

void
add_results_to_results (Results * res, Results * to_add)
{
    res->mobility += to_add->mobility;
    res->diffusivity += to_add->diffusivity;
    res->currentDensity += to_add->currentDensity;
    res->simulationTime += to_add->simulationTime;
    res->totalSimulationTime += to_add->totalSimulationTime;
    res->equilibrationEnergy += to_add->equilibrationEnergy;
    res->avgenergy += to_add->avgenergy;
    res->einsteinrelation += to_add->einsteinrelation;

    res->analytic_mobility += to_add->analytic_mobility;
    res->analytic_fermienergy += to_add->analytic_fermienergy;
    res->analytic_transportenergy += to_add->analytic_transportenergy;

}

void
divide_results_by_scalar (Results * res, double scalar)
{

    res->mobility /= scalar;
    res->diffusivity /= scalar;
    res->currentDensity /= scalar;
    res->simulationTime /= scalar;
    res->totalSimulationTime /= scalar;
    res->equilibrationEnergy /= scalar;
    res->avgenergy /= scalar;
    res->einsteinrelation /= scalar;

    res->analytic_mobility /= scalar;
    res->analytic_fermienergy /= scalar;
    res->analytic_transportenergy /= scalar;

}
