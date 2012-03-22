//  Created by Jan Oliver Oelerich

#include "hop.h"

void al_calcMobility (Results * res);
void al_calcTransportEnergy (Results * res);
void al_calcFermiEnergy (Results * res);
void al_normalizeDos (Results * params);
double al_DOS (double x, Results * params);
double al_DOSunnormalized (double x, void *p);
double al_dosfermi (double x, void *p);
double al_fermidirac (double x, Results * params);
double al_RadiusIntegral (double x, void *p);

// the normalization of the dos function
double dosnorm = 0.0;

// the integration workspace
gsl_integration_workspace *w;

/*
 * This function calculates various properties of the
 * system in the right order. 
 */
void
AL_run (Results * total, int *iRun)
{
    w = gsl_integration_workspace_alloc (1000000);

    al_normalizeDos (&(total[*iRun - 1]));
    al_calcFermiEnergy (&(total[*iRun - 1]));
    al_calcTransportEnergy (&(total[*iRun - 1]));
    al_calcMobility (&(total[*iRun - 1]));

}

/*
 * This function calculates the mobility for the system.
 */
void
al_calcMobility (Results * res)
{

    double error, int1, int2, radius, prefactor, exponent, end;

    end = res->analytic_transportenergy / prms.temperature;

    gsl_function F;

    // integral 1: Radius
    F.function = &al_RadiusIntegral;
    F.params = res;
    gsl_integration_qagiu (&F, end, 0, 1e-7, 100000, w, &int1, &error);

    // integral 2: concentration
    F.function = &al_dosfermi;
    gsl_integration_qagiu (&F, end, 0, 1e-7, 100000, w, &int2, &error);

    // calculate some variables
    prefactor = 4 * M_PI / (3 * prms.percthresh);
    radius = 1 / (pow (prefactor * int1, 1.0 / 3.0));
    exponent =
        -2 * radius / prms.loclength + end -
        res->analytic_fermienergy / prms.temperature;

    // mobility
    res->analytic_mobility =
        exp (exponent - log (prefactor * radius * int2)) / prms.temperature;
}

/*
 * Find the Fermi energy according to the DOS and the
 * filling of the system.
 */
void
al_calcFermiEnergy (Results * res)
{

    double result, error, precision = 0.00000001;
    double upper = 60, lower = 0;

    gsl_function F;
    F.function = &al_dosfermi;
    F.params = res;

    while (upper - lower > precision)
    {
        res->analytic_fermienergy = (upper + lower) / 2.0;
        gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000000, w, &result, &error);

        if (result < (float) prms.ncarriers / prms.nsites)
            upper = res->analytic_fermienergy;
        else
            lower = res->analytic_fermienergy;
    }

    res->analytic_fermienergy *= prms.temperature;
}

/*
 * Normalize the dos
 */
void
al_normalizeDos (Results * res)
{
    double error;

    gsl_function F;
    F.function = &al_DOSunnormalized;

    gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000000, w, &(dosnorm), &error);

}

/*
 * Find the transport energy of the system iteratively.
 */
void
al_calcTransportEnergy (Results * res)
{

    double precision = 0.0000000001;
    double upper = 40.0, lower = 0.0, x;
    double result, error, rightside;

    gsl_function F;
    F.function = &al_RadiusIntegral;
    F.params = res;

    double factor = pow (4 * M_PI / (3 * prms.percthresh),
                         1.0 / 3.0) * prms.loclength;

    // now the right side...
    while (upper - lower > precision)
    {
        x = (upper + lower) / 2.0;

        gsl_integration_qagiu (&F, x, 0, 1e-7, 1000000, w, &result, &error);
        rightside =
            al_RadiusIntegral (x, res) / (pow (result, 4. / 3.) * 1.5 * factor);

        if (rightside > 1.0)
            upper = x;
        else
            lower = x;

    }

    gsl_integration_qagiu (&F, x, 0, 1e-7, 1000000, w, &result, &error);
    rightside =
        al_RadiusIntegral (x, res) / (pow (result, 4. / 3.) * 1.5 * factor);

    if (fabs (rightside - 1.0) > 0.01)
        x = 0;

    res->analytic_transportenergy = x * prms.temperature;
}

// the following are helper functions --------------

double
al_DOSunnormalized (double x, void *p)
{
    return exp (-pow (x / prms.sigma, prms.exponent));
}

double
al_DOS (double x, Results * params)
{
    return al_DOSunnormalized (x, params) / dosnorm;
}

double
al_fermidirac (double x, Results * params)
{
    return 1 / (1 + exp (params->analytic_fermienergy / prms.temperature - x));
}

double
al_dosfermi (double x, void *p)
{
    Results *params = (Results *) p;
    return al_DOS (x, params) * al_fermidirac (x, params);
}

double
al_RadiusIntegral (double x, void *p)
{
    Results *params = (Results *) p;
    return al_DOS (x, params) * (1 - al_fermidirac (x, params));
}
