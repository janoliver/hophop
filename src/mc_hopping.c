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

void hoppingStep (Carrier * carriers, RunParams * runprms);
void hop (Carrier * c, SLE * dest, Vector * dist, RunParams * runprms);
void updateCarrier (Carrier * carriers, RunParams * runprms);

/*
 * This function runs the iteration of the simulation.  It keeps track of
 * the passed "time" (in arbitrary units) and ouputs the progress.
 */
void
MC_simulation (Site * sites, Carrier * carriers, RunParams * runprms,
               int iReRun)
{
    int j;
    int ncarriers = 1;

    // is this the meanfield mode?
    if (prms.many)
        ncarriers = prms.ncarriers;

    // we need the current simulation time.
    double simTimeOld = runprms->simulationTime;
    runprms->nHops = 0;
    runprms->stat = false;

    // relaxation, no time or hop counting
    for (j = 0; j <= 100; j++)
    {
        while (runprms->nHops < prms.relaxation / 100 * j)
            hoppingStep (carriers, runprms);

        output (O_SERIAL, "\r\tRelaxing...   (run %d of %d):\t%2d%%", iReRun,
                prms.number_reruns, (int) j);
        fflush (stdout);
    }
    output (O_SERIAL, " Done.\n");

    // actual simulation, time and hop counting
    // we need to renormalize the carrier occupation time and simulation time, since
    // we might have multiple runs due to -x, while simulationTime is counted totally
    runprms->nHops = 0;
    runprms->stat = true;

    for (j = 0; j < ncarriers; ++j)
        carriers[j].occTime -= (runprms->simulationTime - simTimeOld);
    runprms->simulationTime = simTimeOld;

    for (j = 0; j <= 100; j++)
    {
        while (runprms->nHops <= prms.simulation / 100 * j)
            hoppingStep (carriers, runprms);

        output (O_SERIAL, "\r\tSimulating... (run %d of %d):\t%2d%%", iReRun,
                prms.number_reruns, (int) j);
        fflush (stdout);
    }
    output (O_SERIAL, " Done.\n");

    // finish statistics
    for (j = 0; j < runprms->nSites; ++j)
        if (sites[j].tempOccTime > 0)
            sites[j].totalOccTime +=
                runprms->simulationTime - sites[j].tempOccTime;

    for (j = 0; j < ncarriers; ++j)
    {
        carriers[j].dx2 += pow (carriers[j].ddx, 2.0);
        carriers[j].dy2 += pow (carriers[j].ddy, 2.0);
        carriers[j].dz2 += pow (carriers[j].ddz, 2.0);
    }

}

/*
 * This function executes one step of the simulation. It finds the carrier
 * or occupied site with the lowest occupation time and defines this as
 * the next electron to hop. It then reduces all the other times of all
 * other electrons and finds the destination for the hop. Then, the hop()
 * function is called with these parameters.
 */
void
hoppingStep (Carrier * carriers, RunParams * runprms)
{
    Carrier *c = NULL;
    SLE *dest = NULL;
    double randomHopProb, probSum;
    int i;

    // heap
    c = &carriers[0];

    // determine the next destination site
    randomHopProb = (float) gsl_rng_uniform (runprms->r) * c->site->rateSum;
    probSum = 0.0;
    for (i = 0; ((i < c->site->nNeighbors) && (probSum <= randomHopProb)); ++i)
    {
        dest = &(c->site->neighbors[i]);
        probSum += dest->rate;
    }

    runprms->simulationTime = c->occTime;

    if (dest->s->carrier == NULL)
    {
        // do the hopping and write some statistics
        runprms->nHops++;

        hop (c, dest, &dest->dist, runprms);
    }
    else
    {
        if (runprms->stat)
        {
            runprms->nFailedAttempts++;
            c->nFailedAttempts++;
        }
    }
    updateCarrier (carriers, runprms);
}

/*
 * In this function, one hop is executed. Here, the track-keeping of
 * energies and so on should happen. It updates both sites that take part
 * in the hopping process and all of the sites around these two sites.
 */
void
hop (Carrier * c, SLE * dest, Vector * dist, RunParams * runprms)
{
    Site *orig = c->site;

    // update origin site
    orig->carrier = NULL;

    if (runprms->stat)
    {
        dest->nTransitions++;

        orig->totalOccTime += runprms->simulationTime - orig->tempOccTime;
        orig->tempOccTime = 0.0;
        dest->s->tempOccTime = runprms->simulationTime;

        c->dx += dist->x;
        c->dy += dist->y;
        c->dz += dist->z;

        c->ddx += dist->x;
        c->ddy += dist->y;
        c->ddz += dist->z;
    }

    // update carrier
    c->site = dest->s;

    // update destination site
    dest->s->carrier = c;
    if (runprms->stat && orig->energy < dest->s->energy)
        dest->s->visitedUpward++;
    else if (runprms->stat)
        dest->s->visited++;

}

/*
 * Subtract the `struct timeval' values X and Y, storing the result in
 * RESULT.  Return 1 if the difference is negative, otherwise 0.
 */
int
timeval_subtract (struct timeval *result, struct timeval *y, struct timeval *x)
{
    /*
     * Perform the carry for the later subtraction by updating y. 
     */
    if (x->tv_usec < y->tv_usec)
    {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000)
    {
        int nsec = (x->tv_usec - y->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    /*
     * Compute the time remaining to wait.
     * tv_usec is certainly positive. 
     */
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;

    /*
     * Return 1 if result is negative. 
     */
    return x->tv_sec < y->tv_sec;
}

/*
 * This is a binary tree search for the carrier with the lowest occupation
 * time. It moves it to the top of the carriers list. The last jumped
 * carrier is assigned a new occupation time.
 */
void
updateCarrier (Carrier * c, RunParams * runprms)
{
    int smallest, i = 0;

    c[0].occTime +=
        (float) gsl_ran_exponential (runprms->r, 1.0) / c[0].site->rateSum;

    // children always sit at c[2i+1] and c[2i+2]
    // do this only when the non-meanfield mc simulation is selected
    Carrier tmp;
    while (1 && prms.many)
    {
        smallest = (2 * i + 1 < prms.ncarriers &&
                    c[i].occTime > c[2 * i + 1].occTime) ? 2 * i + 1 : i;
        if (2 * i + 2 < prms.ncarriers &&
            c[smallest].occTime > c[2 * i + 2].occTime)
            smallest = 2 * i + 2;

        if (i == smallest)
            break;

        // swap
        tmp = c[i];
        c[i] = c[smallest];
        c[smallest] = tmp;
        i = smallest;
    }
}
