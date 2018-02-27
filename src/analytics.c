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


typedef struct a_params
{
    double dosnormalization;
    double fermienergy;
} a_params;

double DOSunnormalized (double x, void *p);
double DOS (double x, a_params * params);
double fermidirac (double x, a_params * params);
double dosfermi (double x, void *p);


double
calcFermiEnergy ()
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000000);
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
    F.params = &params;

    while (upper - lower > precision)
    {
        params.fermienergy = (upper + lower) / 2.0;
        if (!prms.gaussian)
            gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000000, w, &result, &error);
        else
            gsl_integration_qagi (&F, 0, 1e-7, 1000000, w, &result, &error);
        if (result < prms.ncarriers * 1.0 / prms.nsites)
            upper = params.fermienergy;
        else
            lower = params.fermienergy;
    }

    return -1.0 * params.fermienergy;
}


double
DOSunnormalized (double x, void *p)
{
    return exp (-pow (x, prms.exponent));
}

double
DOS (double x, a_params * params)
{
    return DOSunnormalized (x, params) / params->dosnormalization;
}

double
fermidirac (double x, a_params * params)
{
    return 1 / (1 + exp ((params->fermienergy - x) / prms.temperature));
}

double
dosfermi (double x, void *p)
{
    struct a_params *params = (struct a_params *) p;
    return DOS (x, params) * fermidirac (x, params);
}
