//  Created by Jan Oliver Oelerich

#include "hop.h"

void
MC_calculateResults (Site * sites, Carrier * carriers, Results * res)
{
    // calculate results
    res->mobility = calcMobility (carriers, res);
    res->diffusivity = calcDiffusivity (carriers, res);
    res->currentDensity = calcCurrentDensity (carriers, res);
    res->equilibrationEnergy = calcEquilibrationEnergy (sites, res);
    res->avgenergy = calcAverageEnergy (carriers);
    res->einsteinrelation = calcEinsteinRelation (carriers, res);
}

/*
 * Here, the carrier mobility in the direction of the electric field is
 * calculated.  Equation: \mu = \frac{<z>}{F*t} where F is the electric
 * field, t is the simulation time and <z> the mean distance in
 * z-direction the elctrons hopped.
 */
double
calcMobility (Carrier * carriers, Results * results)
{
    int i;
    double ez = 0;
    for (i = 0; i < prms.ncarriers; ++i)
        ez += carriers[i].dz;
    return ez / (prms.ncarriers * prms.field * results->simulationTime);
}

/*
 * This is the average energy weighted by occupation time of this energy.
 * It should be something like the equlibration energy.
 */
double
calcEquilibrationEnergy (Site * sites, Results * results)
{
    size_t i;
    double sum = 0;

    for (i = 0; i < prms.nsites; ++i)
        sum += sites[i].totalOccTime * sites[i].energy;

    return sum / (prms.ncarriers * results->simulationTime);
}

/*
 * Here, the diffusivity of the system is calculated.  This is, right now,
 * perpendicular to the field direction.
 */
double
calcDiffusivity (Carrier * carriers, Results * results)
{
    double ex, ey, ez, ex2, ey2, ez2, /*ds, */ dp /*, sez */ ;
    int i;
    ex = 0.0;
    ey = 0.0;
    ez = 0.0;
    ex2 = 0.0;
    ey2 = 0.0;
    ez2 = 0.0;

    for (i = 0; i < prms.ncarriers; ++i)
    {
        ex += carriers[i].dx;
        ey += carriers[i].dy;
        ez += carriers[i].dz;
        ex2 += carriers[i].dx2;
        ey2 += carriers[i].dy2;
        ez2 += carriers[i].dz2;
    }
    ex2 = ex2 / prms.ncarriers;
    ey2 = ey2 / prms.ncarriers;
    ez2 = ez2 / prms.ncarriers;
    ex = ex / prms.ncarriers;
    ey = ey / prms.ncarriers;
    ez = ez / prms.ncarriers;

    // sez = pow(ez / prms.ncarriers, 2.0);
    //ds = (ez2 - sez) / (2 * results->simulationTime);
    dp = (ex2 + ey2) / (4);
    return dp;
}

/*
 * Calculate the einstein relation in units of e/sigma.
 */
double
calcEinsteinRelation (Carrier * carriers, Results * results)
{
    double ex2, ey2, ez;
    int i;
    ez = 0.0;
    ex2 = 0.0;
    ey2 = 0.0;

    for (i = 0; i < prms.ncarriers; ++i)
    {
        ex2 += carriers[i].dx2;
        ey2 += carriers[i].dy2;
        ez += carriers[i].dz;
    }
    ex2 /= prms.ncarriers;
    ey2 /= prms.ncarriers;
    ez /= prms.ncarriers;

    return (4. * ez) / (prms.field * (ex2 + ey2) * results->simulationTime);

}

/*
 * This function determines the current density of the system.
 */
double
calcCurrentDensity (Carrier * carriers, Results * results)
{
    int i;
    double ez = 0;
    for (i = 0; i < prms.ncarriers; ++i)
        ez += carriers[i].dz;
    return ez / (results->simulationTime * prms.nsites);
}

/*
 * The average energy of the charge carriers at the present moment
 * is calculated here. 
 */
double
calcAverageEnergy (Carrier * carriers)
{
    int i;
    double avg = 0;
    for (i = 0; i < prms.ncarriers; ++i)
        avg += carriers[i].site->energy;

    return avg / prms.ncarriers;
}
