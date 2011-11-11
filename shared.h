//  Created by Jan Oliver Oelerich

#ifndef _SHARED_H_
#define _SHARED_H_

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include "cmdline.h"

typedef struct results {
    double mobility;
    double diffusivity;
    double simulationTime;
    double fermiEnergy;
    double transportEnergy;
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




#endif /* _SHARED_H_ */
