//  Created by Jan Oliver Oelerich

#include "hop.h"

// we use this struct to create a temporary linked list
// of the found neighbors.
typedef struct linked_neighbors
{
    Site *s;
    Vector dist;
    struct linked_neighbors *next;
} ln;

// this struct is used to categorize sites
// in space for easier finding of the neighbors
typedef struct cell
{
    int counter;
    int index;
    ln *siteList;
} Cell;

// we use this struct to store softpairs
typedef struct softpair
{
    Site *i;
    Site *j;
    struct softpair *next;
} Softpair;



Cell *createCells (Site * sites);
Cell *getCell3D (Cell * cells, ssize_t x, ssize_t y, ssize_t z);
void setNeighbors (Site * s, Cell * cells);
double calcHoppingRate (Site i, Site j);
Vector distance (Site * i, Site * j);
int compare_neighbors (const void *a, const void *b);

/*
 * creates n sites randomly distributed within a box of the size X*Y*Z
 * cells (!) and with random energies picked out of a gaussian
 * distribution with standard deviation sigma (=1.0 since this should be
 * the energy scale!).
 */
Site *
MC_createSites (RunParams * runprms)
{
    size_t i, j, k, l;
    Site *s;

    s = (Site *) malloc (prms.nsites * sizeof (Site));

    // lattice case. Map site index to x,y,z coordinates
    for (l = 0; l < prms.length_x; ++l)
        for (j = 0; j < prms.length_y; ++j)
            for (k = 0; k < prms.length_z; ++k)
            {
                i = ((l * prms.length_y + j) * prms.length_z) + k;
                if (prms.lattice)
                {
                    s[i].x = (float) l;
                    s[i].y = (float) j;
                    s[i].z = (float) k;
                }
                else
                {
                    s[i].x =
                        (float) gsl_rng_uniform (runprms->r) * prms.length_x;
                    s[i].y =
                        (float) gsl_rng_uniform (runprms->r) * prms.length_y;
                    s[i].z =
                        (float) gsl_rng_uniform (runprms->r) * prms.length_z;
                }

                if (!prms.gaussian)
                    s[i].energy = (float) gsl_ran_exppow (runprms->r, 1.,
                                                          prms.exponent);
                else
                    s[i].energy = (float) gsl_ran_gaussian (runprms->r, 1.);

                s[i].carrier = NULL;
                s[i].visited = 0;
                s[i].visitedUpward = 0;
                s[i].index = i;
                s[i].totalOccTime = 0.0;
                s[i].tempOccTime = 0.0;
                s[i].neighbors = NULL;
                s[i].nNeighbors = 0;
                s[i].rateSum = 0.0;
            }

    return s;
}

/*
 * allocate carrier array and initialize the values
 */
Carrier *
MC_createCarriers (Site * sites)
{
    int i;
    Carrier *c;

    // allocate carrier memory
    c = (Carrier *) malloc (sizeof (Carrier) * prms.ncarriers);

    // distribute carriers 
    for (i = 0; i < prms.ncarriers; ++i)
    {
        c[i].index = i;
        c[i].dx = 0.0;
        c[i].dy = 0.0;
        c[i].dz = 0.0;
        c[i].dx2 = 0.0;
        c[i].dy2 = 0.0;
        c[i].dz2 = 0.0;
        c[i].nFailedAttempts = 0;
        c[i].nHops = 0;
    }

    return c;
}


/*
 * Distributes the electrons over random sites. Since all site energies
 * and positions were created randomly, simply the first nCarriers sites
 * can be set to occupied.m This can be modified to have a localized
 * source of electrons or anything else.
 */
void
MC_distributeCarriers (Carrier * c, Site * sites, RunParams * runprms,
                       Results * res)
{
    size_t i;
    Carrier tmp;

    // select prms.ncarriers random sites for the carriers
    Site *sample = malloc (prms.ncarriers * sizeof (Site));
    gsl_ran_choose (runprms->r, sample, prms.ncarriers, sites, prms.nsites,
                    sizeof (Site));

    // clear all the sites and set the correct index
    for (i = 0; i < prms.nsites; ++i)
    {
        sites[i].carrier = NULL;
        sites[i].tempOccTime = 0;
        sites[i].index = i;
    }

    // distribute carriers 
    for (i = 0; i < prms.ncarriers; ++i)
    {
        c[i].site = &sites[sample[i].index];
        c[i].occTime = res->simulationTime +
            (float) gsl_ran_exponential (runprms->r, 1.0) / c[i].site->rateSum;

        // now insert the carrier into the heap
        while (i > 0 && c[i].occTime < c[i / 2].occTime)
        {
            tmp = c[i];
            c[i] = c[i / 2];
            c[i / 2] = tmp;
            i /= 2;
        }
        c[i].ddx = 0.0;
        c[i].ddy = 0.0;
        c[i].ddz = 0.0;
        c[i].site->carrier = &c[i];
        c[i].site->tempOccTime = 0.000001;
    }
}

/*
 * Divides the sample into equally sized cells, so that it's easier
 * to keep track of relevant neighbors. One Cell has the size
 * 2 * rc * a, where rc is the cut-off radius for the hops, 
 * a is the localization length of the WF. 
 * It also stores all the sites that lay within the cell
 * into a linked list which is a member of the cell struct.
 */
Cell *
createCells (Site * sites)
{
    size_t i, tx, ty, tz, nCells;
    Cell *c, *temp;
    ln *head;

    nCells = prms.nx * prms.ny * prms.nz;

    c = (Cell *) malloc (nCells * sizeof (Cell));

    for (i = 0; i < nCells; ++i)
    {
        c[i].counter = 0;
        c[i].siteList = NULL;
        c[i].index = i;
    }

    for (i = 0; i < prms.nsites; ++i)
    {
        tx = sites[i].x / prms.cutoff_radius;
        ty = sites[i].y / prms.cutoff_radius;
        tz = sites[i].z / prms.cutoff_radius;

        temp = getCell3D (c, tx, ty, tz);
        if (temp->counter == 0)
            temp->siteList = NULL;

        temp->counter++;

        head = temp->siteList;
        temp->siteList = (ln *) malloc (sizeof (ln));
        temp->siteList->s = &sites[i];
        temp->siteList->next = head;
    }

    return c;
}

/*
 * This function creates the array Site.neighbors, assigning neighbors of
 * the sites and the corresponding hopping rate to it. It takes a while,
 * especially for big values of nSites.  It is also the main consumer of
 * memory!!
 */
void
MC_createHoppingRates (Site * sites)
{
    size_t k, l;
    ln *sList, *tmp;
    Cell *cells = (Cell *) createCells (sites);

    for (l = 0; l < 100; ++l)
    {
        // this weirdness with 99 takes care of site numbers that
        // cannot be divided by 100.
        for (k = l * prms.nsites / 99; k < ((l + 1) * prms.nsites / 99); ++k)
            if (k < prms.nsites)
                setNeighbors (&sites[k], cells);

        output (O_SERIAL, "\r\tInitializing...: \t\t%2d%%", (int) l);
        fflush (stdout);
    }
    output (O_SERIAL, " Done.\n");

    // free cells
    for (l = 0; l < prms.nx * prms.ny * prms.nz; ++l)
    {
        sList = cells[l].siteList;
        while (sList)
        {
            tmp = sList->next;
            free (sList);
            sList = tmp;
        }
    }

    free (cells);
}

/*
 * This function removes the found softpairs, using the algorithm
 * described in the PhD Thesis of Fredrik Jansson based ok xyz.
 */
void
MC_removeSoftPairs (Site * sites)
{
    Softpair *softpair, *temp, *newSoftpair = NULL, *sp = NULL;
    int j, i, nSoftPairs = 0;
    bool ignore;
    SLE *ntemp = NULL, *neighbor = NULL;
    float rateb, rateab, ratebx = 0.0, rateSum;

    // find softpairs
    for (j = 0; j < prms.nsites; ++j)
    {
        for (i = 0; i < sites[j].nNeighbors; ++i)
        {
            if (sites[j].neighbors[i].rate / sites[j].rateSum >
                prms.softpairthreshold)
            {
                newSoftpair = (Softpair *) malloc (sizeof (Softpair));
                newSoftpair->i = &sites[j];
                newSoftpair->j = sites[j].neighbors[i].s;
                newSoftpair->next = sp;
                sp = newSoftpair;
                nSoftPairs++;
            }
        }
    }

    if (nSoftPairs > 0)
    {
        // some output
        output (O_SERIAL, "\tRemoving Softpairs...:");

        temp = newSoftpair;
        softpair = newSoftpair;

        // remove softpairs
        while (softpair)
        {
            // check if pair has been found twice. If so, skip it.
            ignore = false;
            while (temp)
            {
                if (softpair->i->index == temp->j->index &&
                    softpair->j->index == temp->i->index &&
                    softpair->j->index > temp->i->index)
                    ignore = true;
                temp = temp->next;
            }

            // remove the softpair
            if (!ignore)
            {
                rateb = 0;
                rateab = 0;

                // find the partner and elimintate the transition
                for (i = 0; i < softpair->i->nNeighbors; ++i)
                {
                    if (softpair->j->index ==
                        softpair->i->neighbors[i].s->index)
                    {
                        neighbor = &(softpair->i->neighbors[i]);

                        rateab = neighbor->rate;
                        rateb = neighbor->s->rateSum;
                        neighbor->rate = 0.0;
                        break;
                    }
                }

                // update transition rates of all other neighbors
                rateSum = 0.0;
                for (i = 0; i < softpair->i->nNeighbors; ++i)
                {
                    neighbor = &(softpair->i->neighbors[i]);

                    // find the partner
                    if (softpair->j->index != neighbor->s->index)
                    {
                        ratebx = 0.0;
                        for (j = 0; j < softpair->j->nNeighbors; ++j)
                        {
                            if (softpair->j->neighbors[j].s->index ==
                                neighbor->s->index)
                                ratebx = ntemp->rate;
                        }
                        if (ratebx > 0)
                        {
                            neighbor->rate *= rateb;
                            neighbor->rate += ratebx * rateab;
                            neighbor->rate /= (rateb + rateab);
                        }
                        rateSum += neighbor->rate;
                    }
                }
                softpair->i->rateSum = rateSum;
            }
            softpair = softpair->next;
        }

        // free memory
        softpair = newSoftpair;
        while (softpair)
        {
            temp = softpair->next;
            free (softpair);
            softpair = temp;
        }


        // all removed ?
        int counter = 0;
        for (j = 0; j < prms.nsites; ++j)
            for (i = 0; i < sites[j].nNeighbors; ++i)
                if (sites[j].neighbors[i].rate / sites[j].rateSum >
                    prms.softpairthreshold)
                    counter++;

        // some output
        output (O_SERIAL, "\t\tDone. %d Softpairs found.\n", counter);

        // recursive call if there are still softpairs
        if (counter > 0)
            MC_removeSoftPairs (sites);
    }
}

/*
 * Returns a linked list of type SLE for the neighbors of the site
 * s. Requires the array of cells.
 */
void
setNeighbors (Site * s, Cell * cells)
{
    int i, k, l;
    Cell *c;
    ln *siteList;
    Vector d;
    double length;

    ln *curr, *head = NULL, *tmp;

    // now, find all the neighbors
    // iterate over all neighboring cells
    for (i = -1; i <= 1; i++)
        for (k = -1; k <= 1; k++)
            for (l = -1; l <= 1; l++)
            {

                c = getCell3D (cells, i + floor (s->x / prms.cutoff_radius),
                               floor (k + s->y / prms.cutoff_radius),
                               floor (l + s->z / prms.cutoff_radius));
                siteList = c->siteList;

                while (siteList)
                {
                    // if this is not the same Site as s and
                    // is within the cutoff radius, attach it 
                    // to the neighbor list and count one up.
                    d = distance (s, siteList->s);
                    length =
                        sqrt (pow (d.x, 2.) + pow (d.y, 2.) + pow (d.z, 2.));

                    if (s->index != siteList->s->index &&
                        length < prms.cutoff_radius)
                    {
                        curr = (ln *) malloc (sizeof (ln));

                        s->nNeighbors++;

                        curr->dist = d;
                        curr->s = siteList->s;
                        curr->next = head;
                        head = curr;

                    }
                    siteList = siteList->next;
                }
            }

    // now make the neighbor array
    curr = head;
    s->neighbors = (SLE *) malloc (sizeof (SLE) * s->nNeighbors);
    i = 0;
    while (curr)
    {
        s->neighbors[i].s = curr->s;
        s->neighbors[i].rate = calcHoppingRate (*s, *curr->s);
        s->neighbors[i].dist = curr->dist;
        s->neighbors[i].nTransitions = 0;
        s->rateSum += s->neighbors[i].rate;

        curr = curr->next;
        ++i;
    }

    // free the temporary list
    curr = head;
    while (curr)
    {
        tmp = curr->next;
        free (curr);
        curr = tmp;
    }

    // sort the neighbors according to the rate to save computation
    // time while simulating
    qsort (s->neighbors, s->nNeighbors, sizeof (SLE), compare_neighbors);
}

/*
 * This function calculates the hopping rate from one site to another. It
 * takes into account periodic boundary con- ditions. float instead of
 * double to save memory.
 */
double
calcHoppingRate (Site i, Site j)
{
    double r = 1.0;
    float dE, dist;
    Vector distances;

    // calc spatial and energetic distances
    distances = distance (&i, &j);
    dE = j.energy - i.energy - prms.field * distances.z;
    dist = sqrt (pow (distances.x, 2.0) +
                 pow (distances.y, 2.0) + pow (distances.z, 2.0));

    // spatial part
    r *= exp (-2.0 * dist / prms.loclength);

    // energy part
    if (dE > 0 && prms.temperature > 0)
        r *= exp (-1.0 * dE / prms.temperature);
    if (dE > 0 && prms.temperature == 0)
        r = 0;

    return r;
}


/*
 * This function calculates the distance between two sites and takes into
 * account periodic boundary conditions. It returns a Vector struct with
 * the fields x,y and z.
 */
Vector
distance (Site * i, Site * j)
{
    float lx, ly, lz;
    Vector vec;

    lx = prms.length_x / 2.;
    ly = prms.length_y / 2.;
    lz = prms.length_z / 2.;

    vec.x = j->x - i->x;
    vec.y = j->y - i->y;
    vec.z = j->z - i->z;

    if (vec.x > lx)
        vec.x -= prms.length_x;
    else if (vec.x < -lx)
        vec.x += prms.length_x;

    if (vec.y > ly)
        vec.y -= prms.length_y;
    else if (vec.y < -ly)
        vec.y += prms.length_y;

    if (vec.z > lz)
        vec.z -= prms.length_z;
    else if (vec.z < -lz)
        vec.z += prms.length_z;

    return vec;
}


/*
 * This function maps a 3D matrix for the cells of the sample to a 1D
 * array (or the other way round.). Expects coordinates of the cell (small
 * x,y,z) as well as the cell array.  Returns a pointer to the desired
 * cell element.
 */
Cell *
getCell3D (Cell * cells, ssize_t x, ssize_t y, ssize_t z)
{
    while (x >= prms.nx)
        x -= prms.nx;
    while (y >= prms.ny)
        y -= prms.ny;
    while (z >= prms.nz)
        z -= prms.nz;
    while (x < 0)
        x += prms.nx;
    while (y < 0)
        y += prms.ny;
    while (z < 0)
        z += prms.nz;

    return &cells[(x * prms.nx + y) * prms.nz + z];
}

/*
 * Compare two neighbors by their rates
 */
int
compare_neighbors (const void *a, const void *b)
{
    double diff = (((SLE *) a)->rate - ((SLE *) b)->rate);
    return diff < 0 ? 1 : (diff > 0) ? -1 : 0;
}
