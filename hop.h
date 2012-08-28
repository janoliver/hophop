//  Created by Jan Oliver Oelerich

#ifndef HOP_H
#define HOP_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdarg.h>

#include <sys/time.h>

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

#include "config.h"

#include "cmdline.h"

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
    float add_to_energy;
    float add_to_number;
    bool addto_dos;

    // balance eq.
    bool balance_eq;
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
} RunParams;

struct site_list_element;
struct site;

typedef struct carrier
{
    double dx, dy, dz, dx2, dy2, dz2, ddx, ddy, ddz;
    struct site *site;
    int index;
    double occTime;
    size_t nFailedAttempts;
    size_t nHops;
} Carrier;

typedef struct site
{
    float energy;
    float x;
    float y;
    float z;
    int visited;
    int visitedUpward;
    Carrier *carrier;
    int index;
    struct site_list_element *neighbors;
    int nNeighbors;
    double rateSum;
    float totalOccTime;
    float tempOccTime;
} Site;

typedef struct results
{
    double mobility;
    double diffusivity;
    double einsteinrelation;
    double simulationTime;
    double currentDensity;
    double equilibrationEnergy;
    size_t nHops;
    size_t nFailedAttempts;
    double avgenergy;
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
void average_results (Results * res, Results * total, Results * error);
void init_results (Results * res);
void add_results_to_results (Results * res, Results * to_add);
void divide_results_by_scalar (Results * res, double scalar);
int output (int mode, const char *fmt, ...);

// output.c
void writeSites (Site * sites, int iRun);
void writeSitesConfig (Site * sites, int iRun);
void writeTransitions (Site * sites, int iRun);
void writeConfig (int iRun);
void writeResults (Results * res, int iRun);
void writeSummary (Results * res, Results * error);

// params
bool strArgGiven (char *arg);
void generateParams (Params * prms, int argc, char **argv);

//mc
void MC_simulation (Site * sites, Carrier * carriers, Results * res,
                    RunParams * runprms, int *iRun, int iReRun);
Site *MC_createSites (RunParams * runprms);
void MC_distributeCarriers (Carrier * carriers, Site * sites,
                            RunParams * runprms, Results * res);
Carrier *MC_createCarriers (Site * sites);
void MC_createHoppingRates (Site * sites);
void MC_removeSoftPairs (Site * sites);
void MC_calculateResults (Site * sites, Carrier * carriers, Results * res);
void MC_run (Results * total, RunParams * runprms, int *iRun);

int timeval_subtract (struct timeval *result,
                      struct timeval *x, struct timeval *y);
double calcMobility (Carrier * carriers, Results * results);
double calcDiffusivity (Carrier * carriers, Results * results);
double calcEinsteinRelation (Carrier * carriers, Results * results);
double calcCurrentDensity (Carrier * carriers, Results * results);
double calcEquilibrationEnergy (Site * sites, Results * results);
double calcAverageEnergy (Carrier * carriers);

// balance equations
void BE_run (Results * total, RunParams * runprms, int *iRun);
void BE_solve (Site * sites, Results * res, int *iRun);


#endif /* HOP_H */
