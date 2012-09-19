#include "hop.h"


typedef struct a_params 
{    
    double dosnormalization;
    double fermienergy;
} a_params;

double DOSunnormalized(double x, void * p);
double DOS(double x, a_params * params);
double fermidirac(double x, a_params * params);
double dosfermi(double x, void * p);


double
calcFermiEnergy()
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000000);
    a_params params;

    // find the norm of the DOS
    double error;
    gsl_function F;
    F.function = &DOSunnormalized;
    
    if (!prms.gaussian)
        gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000000, w,
            &(params.dosnormalization), &error);
    else
        gsl_integration_qagi (&F, 0, 1e-7, 1000000, w,
            &(params.dosnormalization), &error);

    double result, precision = 0.00000001;
    double upper = 60, lower = 0;

    F.function = &dosfermi;
    F.params   = &params;

    while(upper-lower > precision)
    {
        params.fermienergy = (upper+lower)/2.0;
        if (!prms.gaussian)
            gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000000, w, &result, &error);
        else
            gsl_integration_qagi (&F, 0, 1e-7, 1000000, w, &result, &error);
        if(result < prms.ncarriers * 1.0 / prms.nsites)
            upper = params.fermienergy;
        else
            lower = params.fermienergy;
    }

    return -1.0 * params.fermienergy;
}


double
DOSunnormalized(double x, void * p)
{
    return exp(-pow(x,prms.exponent));
}

double
DOS(double x, a_params * params)
{
    return DOSunnormalized(x,params)/params->dosnormalization;
}

double 
fermidirac(double x, a_params * params)
{
    return 1 / (1 + exp((params->fermienergy - x)/prms.temperature));
}

double
dosfermi(double x, void * p)
{
    struct a_params * params = (struct a_params *) p;
    return DOS(x, params) * fermidirac(x, params);
}
