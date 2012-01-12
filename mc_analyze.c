//  Created by Jan Oliver Oelerich

#include "mc.h"

void
MC_calculateResults(Site * sites, Carrier * carriers, Results * res)
{
    // calculate results
    res->mobility        = calcMobility(carriers, res);
    res->diffusivity     = calcDiffusivity(carriers, res);
    res->fermiEnergy     = calcFermiEnergy(sites, res);
    res->transportEnergy = calcTransportEnergy(sites, res);
    res->currentDensity  = calcCurrentDensity(carriers, res);
    res->equilibrationEnergy = calcEquilibrationEnergy(sites, res);
}

/**
    Here, the carrier mobility in the direction of the electric field is
    calculated.  Equation: \mu = \frac{<z>}{F*t} where F is the electric
    field, t is the simulation time and <z> the mean distance in
    z-direction the elctrons hopped.
 */
double
calcMobility(Carrier * carriers, Results * results)
{
    int i;
    double ez = 0;
    for(i = 0; i < args.ncarriers_arg; ++i)
        ez += carriers[i].dz;
    
    return ez / (args.ncarriers_arg *
                 args.field_arg * results->simulationTime);
}

/**
   This is the average energy weighted by occupation time of this energy.
   It should be something like the equlibration energy.
 */
double
calcEquilibrationEnergy(Site * sites, Results * results)
{
    size_t i;
    double sum = 0;
    
    for(i = 0; i < args.nsites_arg; ++i)
        sum += sites[i].totalOccTime * sites[i].energy;
   
    return sum / (args.ncarriers_arg * results->simulationTime);
}

/**
    Here, the diffusivity of the system is calculated.  This is, right now,
    perpendicular to the field direction.
 */
double
calcDiffusivity(Carrier * carriers, Results * results)
{
    double ex, ey, ez, ex2, ey2, ez2, /*ds,*/ dp/*, sez*/;
    int i;
    ex  = 0.0;
    ey  = 0.0;
    ez  = 0.0;
    ex2 = 0.0;
    ey2 = 0.0;
    ez2 = 0.0;

    for(i = 0; i < args.ncarriers_arg; ++i)
    {
        ex += carriers[i].dx;
        ey += carriers[i].dy;
        ez += carriers[i].dz;
        ex2 += pow(carriers[i].dx, 2.0);
        ey2 += pow(carriers[i].dy, 2.0);
        ez2 += pow(carriers[i].dz, 2.0);
    }
    ex2 = ex2 / args.ncarriers_arg;
    ey2 = ey2 / args.ncarriers_arg;
    ez2 = ez2 / args.ncarriers_arg;
    ex = ex / args.ncarriers_arg;
    ey = ey / args.ncarriers_arg;
    ez = ez / args.ncarriers_arg;

    // sez = pow(ez / args.ncarriers_arg, 2.0);
    //ds = (ez2 - sez) / (2 * results->simulationTime);
    dp = (ex2 + ey2) / (4 * results->simulationTime);
    return dp;
}

/**
    This function determines the fermi energy of the system.
 */
double
calcFermiEnergy(Site * sites, Results * results)
{
    int i, m;
    double totalEnergy = 0.0;
    double timeRate;
    m = 0;

    for(i = 0; i < args.nsites_arg; ++i)
    {
        timeRate = (sites[i].totalOccTime / results->simulationTime);
        if(0.45 < timeRate && timeRate < 0.52)
        {
            m++;
            totalEnergy += sites[i].energy;
        }
    }
    return totalEnergy / m;
}

/**
    This function determines the transport energy of the system.
 */
double
calcTransportEnergy(Site * sites, Results * results)
{
    return 0.0;
}

/**
    This function determines the current density of the system.
 */
double
calcCurrentDensity(Carrier * carriers, Results * results)
{
    int i;
    double ez = 0;
    for(i = 0; i < args.ncarriers_arg; ++i)
        ez += carriers[i].dz;
    return ez / (results->simulationTime * args.X_arg * args.Y_arg * args.Z_arg);
}
