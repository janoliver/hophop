//  Created by Jan Oliver Oelerich

#ifndef MC_H
#define MC_H

#include "hop.h"


typedef struct vector {
    float x,y,z;
} Vector;

typedef struct site_list_element {
    struct site_list_element * next;
    struct site * s;
    double rate;
    Vector dist;
    int nTransitions;
} SLE;

struct softpair;
typedef struct softpair {
    Site * i;
    Site * j;
    struct softpair * next;
} Softpair;

// hopping_mc.c
int timeval_subtract(struct timeval * result,
                     struct timeval * x, struct timeval * y);

// analyze_mc.c
float calcMobility(Carrier * carriers, Results * results);
float calcDiffusivity(Carrier * carriers, Results * results);
float calcFermiEnergy(Site * sites, Results * results);
float calcTransportEnergy(Site * sites, Results * results);
float calcCurrentDensity(Carrier * carriers, Results * results);


#endif /* MC_H */
