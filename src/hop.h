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


#ifndef HOP_H
#define HOP_H

#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stdlib.h>

#include <sys/time.h>
#include <sys/types.h>

#include <time.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_mksa.h>

#include "cli/cmdline.h"

// the output modes
#define O_PARALLEL 1
#define O_SERIAL   2
#define O_BOTH     3
#define O_FORCE    0

typedef struct params
{
    // all parameters here
    int length_x, length_y, length_z;
    int nx, ny, nz;
    int ncarriers;
    bool many;
    int nsites;
    float exponent;
    float loclength;
    float cutoff_radius;
    float field;
    float temperature;
    bool gaussian;
    bool lattice;
    long relaxation;
    long simulation;
    bool removesoftpairs;
    float softpairthreshold;
    int number_runs;
    int number_reruns;
    bool parallel;
    bool quiet;
    char *output_folder;
    bool output_transitions;
    char *output_summary;
    char *comment;
    bool memreq;
    int nthreads;

    // cutting out and adding sites
    float cut_out_energy;
    float cut_out_width;
    bool cut_dos;

    // balance eq.
    bool balance_eq;
    bool mgmres;
    bool lis;
    float be_abs_tol;
    float be_rel_tol;
    int be_it;
    int be_outer_it;

    // random number stuff
    long rseed;
    time_t curtime;
    struct tm *loctime;
    const gsl_rng_type *T;

    struct gengetopt_args_info *cmdlineargs;

} Params;

// this struct is instantiated for each run of the simulation, also in
// parallel mode. This is done so that the RNG for example is not shared
// between runs. 
typedef struct run_params
{
    gsl_rng *r;
    long rseed_used;

    double simulationTime;
    long nHops;
    long nFailedAttempts;
    int nSites;
    int iRun;
    bool stat;
} RunParams;

struct site_list_element;
struct site;

typedef struct carrier
{
    double dx, dy, dz, dx2, dy2, dz2, ddx, ddy, ddz;
    struct site *site;
    int index;
    double occTime;
    long nFailedAttempts;
    long nHops;
} Carrier;

typedef struct site
{
    float energy;
    float x;
    float y;
    float z;
    unsigned long visited;
    unsigned long visitedUpward;
    Carrier *carrier;
    int index;
    struct site_list_element *neighbors;
    int nNeighbors;
    double rateSum;
    float totalOccTime;
    float tempOccTime;
} Site;

typedef struct result
{
    double *values;
    double avg;
    bool *done;
    double err;
} Result;

typedef struct results_struct
{
    Result mobility;
    Result diffusivity;
    Result einsteinrelation;
    Result currentDensity;
    Result equilibrationEnergy;
    Result avgenergy;

    Result nHops;
    Result simulationTime;
    Result nFailedAttempts;
    Result nSites;

    time_t time_start;
    time_t time_finished;
} Results;


typedef struct vector
{
    float x, y, z;
} Vector;

typedef struct site_list_element
{
    Site *s;
    double rate;
    Vector dist;
    int nTransitions;
} SLE;

extern Params prms;

// helpers
void average_errors (Results * res);
void free_results (Results * res);
void init_results (Results * res);
int output (int mode, const char *fmt, ...);

// output.c
void writeSites (Site * sites, RunParams * runprms);
void writeTransitions (Site * sites, RunParams * runprms);
void writeConfig (RunParams * runprms);
void writeResults (Results * res, RunParams * runprms);
void writeSummary (Results * res);

// params
bool strArgGiven (char *arg);
void generateParams (Params * prms, int argc, char **argv);

//mc
void MC_simulation (Site * sites, Carrier * carriers, RunParams * runprms,
                    int iReRun);
Site *MC_createSites (RunParams * runprms);
void MC_distributeCarriers (Carrier * carriers, Site * sites,
                            RunParams * runprms);
Carrier *MC_createCarriers ();
void MC_createHoppingRates (Site * sites, RunParams * runprms);
void MC_removeSoftPairs (Site * sites, RunParams * runprms);
void MC_calculateResults (Site * sites, Carrier * carriers, Results * res,
                          RunParams * runprms);
void MC_run (Results * total, RunParams * runprms);

int timeval_subtract (struct timeval *result,
                      struct timeval *x, struct timeval *y);

// balance equations
void BE_run (Results * res, RunParams * runprms);

// analytics
double calcFermiEnergy ();

#endif /* HOP_H */
