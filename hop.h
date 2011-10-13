//  Created by Jan Oliver Oelerich

#ifndef HOP_H
#define HOP_H

#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "cmdline.h"

struct site_list_element;
struct site;

typedef struct carrier {
    double dx, dy, dz;
    struct site * site;
    int index;
    double occTime;
    size_t nFailedAttempts;
    size_t nHops;
} Carrier;

typedef struct site {
    float energy;
    float x;
    float y;
    float z;
    int visited;
    int visitedUpward;
    Carrier * carrier;
    int index;
    struct site_list_element * neighbors;
    int nNeighbors;
    double rateSum;
    //double occTime;
    double totalOccTime;
    double tempOccTime;
    double updatedAt;
} Site;

typedef struct results {
    double mobility;
    float  diffusivity;
    double simulationTime;
    float  fermiEnergy;
    float  transportEnergy;
    double currentDensity;
    size_t nHops;
    size_t nFailedAttempts;
} Results;

typedef signed long ssize_t;

extern struct gengetopt_args_info args;
extern unsigned long RSEED;
extern time_t curtime;
struct tm *loctime;
extern const gsl_rng_type * T;
extern gsl_rng * r;

// number of cells in each spatial direction
extern int nx, ny, nz;

void 	  MC_simulation(Site * sites, Carrier * carriers, Results * res);
Site *    MC_createSites();
Carrier * MC_distributeCarriers(Site * sites);
void 	  MC_createHoppingRates(Site * sites);
void      MC_removeSoftPairs(Site * sites);
void 	  MC_calculateResults(Site * sites, Carrier * carriers, Results * res);

// output.c
void writeSites(Site * sites);
void writeSitesConfig(Site * sites);
void writeConfig();
void writeResults(Results * res);
void writeSummary(Results * res, Results * error);

#endif /* HOP_H */
