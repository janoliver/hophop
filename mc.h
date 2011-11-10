//  Created by Jan Oliver Oelerich

#ifndef MC_H
#define MC_H

#include "hop.h"

// hopping_mc.c
int timeval_subtract(struct timeval * result,
                     struct timeval * x, struct timeval * y);

// analyze_mc.c
double calcMobility(Carrier * carriers, Results * results);
double calcDiffusivity(Carrier * carriers, Results * results);
double calcFermiEnergy(Site * sites, Results * results);
double calcTransportEnergy(Site * sites, Results * results);
double calcCurrentDensity(Carrier * carriers, Results * results);


#endif /* MC_H */
