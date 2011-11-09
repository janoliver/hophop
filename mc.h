//  Created by Jan Oliver Oelerich

#ifndef MC_H
#define MC_H

#include "hop.h"

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
