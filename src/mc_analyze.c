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

double calcMobility (Carrier * carriers, RunParams * runprms);
double calcDiffusivity (Carrier * carriers);
double calcEinsteinRelation (Carrier * carriers, RunParams * runprms);
double calcCurrentDensity (Carrier * carriers, RunParams * runprms);
double calcEquilibrationEnergy (Site * sites, RunParams * runprms);
double calcAverageEnergy (Carrier * carriers);

void
MC_calculateResults (Site * sites, Carrier * carriers, Results * res,
                     RunParams * runprms)
{
    // calculate results
    res->mobility.values[runprms->iRun - 1] = calcMobility (carriers, runprms);
    res->mobility.done[runprms->iRun - 1] = true;

    res->diffusivity.values[runprms->iRun - 1] =
        calcDiffusivity (carriers);
    res->diffusivity.done[runprms->iRun - 1] = true;

    res->currentDensity.values[runprms->iRun - 1] =
        calcCurrentDensity (carriers, runprms);
    res->currentDensity.done[runprms->iRun - 1] = true;

    res->equilibrationEnergy.values[runprms->iRun - 1] =
        calcEquilibrationEnergy (sites, runprms);
    res->equilibrationEnergy.done[runprms->iRun - 1] = true;

    res->avgenergy.values[runprms->iRun - 1] = calcAverageEnergy (carriers);
    res->avgenergy.done[runprms->iRun - 1] = true;

    res->einsteinrelation.values[runprms->iRun - 1] =
        calcEinsteinRelation (carriers, runprms);
    res->einsteinrelation.done[runprms->iRun - 1] = true;

    res->simulationTime.values[runprms->iRun - 1] = runprms->simulationTime;
    res->simulationTime.done[runprms->iRun - 1] = true;

    res->nHops.values[runprms->iRun - 1] = runprms->nHops;
    res->nHops.done[runprms->iRun - 1] = true;

    res->nFailedAttempts.values[runprms->iRun - 1] = runprms->nFailedAttempts;
    res->nFailedAttempts.done[runprms->iRun - 1] = true;

    res->nSites.values[runprms->iRun - 1] = (float) runprms->nSites;
    res->nSites.done[runprms->iRun - 1] = true;

}

/*
 * Here, the carrier mobility in the direction of the electric field is
 * calculated.  Equation: \mu = \frac{<z>}{F*t} where F is the electric
 * field, t is the simulation time and <z> the mean distance in
 * z-direction the elctrons hopped.
 */
double
calcMobility (Carrier * carriers, RunParams * runprms)
{
    int i;
    double ez = 0;

    // the meanfield stuff
    int ncarriers = 1;
    if (prms.many)
        ncarriers = prms.ncarriers;

    for (i = 0; i < ncarriers; ++i)
        ez += carriers[i].dz;
    return ez / (ncarriers * prms.field * runprms->simulationTime);
}

/*
 * This is the average energy weighted by occupation time of this energy.
 * It should be something like the equlibration energy.
 */
double
calcEquilibrationEnergy (Site * sites, RunParams * runprms)
{
    long i;
    double sum = 0;

    // the meanfield stuff
    int ncarriers = 1;
    if (prms.many)
        ncarriers = prms.ncarriers;

    for (i = 0; i < runprms->nSites; ++i)
        sum += sites[i].totalOccTime * sites[i].energy;

    return sum / (ncarriers * runprms->simulationTime);
}

/*
 * Here, the diffusivity of the system is calculated.  This is, right now,
 * perpendicular to the field direction.
 */
double
calcDiffusivity (Carrier * carriers)
{
    double ex2, ey2;
    int i;
    ex2 = 0.0;
    ey2 = 0.0;

    // the meanfield stuff
    int ncarriers = 1;
    if (prms.many)
        ncarriers = prms.ncarriers;

    for (i = 0; i < ncarriers; ++i)
    {
        ex2 += carriers[i].dx2;
        ey2 += carriers[i].dy2;
    }
    ex2 /= (ncarriers);
    ey2 /= (ncarriers);
    return (ex2 + ey2) / (4);

}

/*
 * Calculate the einstein relation in units of e/sigma.
 */
double
calcEinsteinRelation (Carrier * carriers, RunParams * runprms)
{
    double ex2, ey2, ez;
    int i;
    ez = 0.0;
    ex2 = 0.0;
    ey2 = 0.0;

    // the meanfield stuff
    int ncarriers = 1;
    if (prms.many)
        ncarriers = prms.ncarriers;

    for (i = 0; i < ncarriers; ++i)
    {
        ex2 += carriers[i].dx2;
        ey2 += carriers[i].dy2;
        ez += carriers[i].dz;
    }

    return (4. * ez) / (runprms->simulationTime * prms.field * (ex2 + ey2));

}

/*
 * This function determines the current density of the system.
 */
double
calcCurrentDensity (Carrier * carriers, RunParams * runprms)
{
    int i;
    double ez = 0;

    // the meanfield stuff
    int ncarriers = 1;
    if (prms.many)
        ncarriers = prms.ncarriers;

    for (i = 0; i < ncarriers; ++i)
        ez += carriers[i].dz;
    return ez / (runprms->simulationTime * runprms->nSites);
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

    // the meanfield stuff
    int ncarriers = 1;
    if (prms.many)
        ncarriers = prms.ncarriers;

    for (i = 0; i < ncarriers; ++i)
        avg += carriers[i].site->energy;

    return avg / ncarriers;
}
