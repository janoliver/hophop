//  Created by Jan Oliver Oelerich

#include "hop.h"

// this struct is used to categorize sites
// in space for easier finding of the neighbors
typedef struct cell
{
    int counter;
    int index;
    SLE *siteList;
} Cell;

Cell *createCells (Site * sites);
Cell *getCell3D (Cell * cells, ssize_t x, ssize_t y, ssize_t z);
void setNeighbors (Site * s, Cell * cells);
double calcHoppingRate (Site i, Site j);
Vector distance (Site * i, Site * j);
SLE *sortNeighbors (SLE * list);

/*
 * creates n sites randomly distributed within a box of the size X*Y*Z
 * cells (!) and with random energies picked out of a gaussian
 * distribution with standard deviation sigma (=1.0 since this should be
 * the energy scale!).
 */
Site *
MC_createSites ()
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
                    s[i].x = (float) gsl_rng_uniform (prms.r) * prms.length_x;
                    s[i].y = (float) gsl_rng_uniform (prms.r) * prms.length_y;
                    s[i].z = (float) gsl_rng_uniform (prms.r) * prms.length_z;
                }

                if (!prms.gaussian)
                    s[i].energy = (float) gsl_ran_exppow (prms.r, prms.sigma,
                                                          prms.exponent);
                else
                    s[i].energy = (float) gsl_ran_gaussian (prms.r, prms.sigma);

                s[i].carrier = NULL;
                s[i].visited = 0;
                s[i].visitedUpward = 0;
                s[i].index = i;
                s[i].totalOccTime = 0.0;
                s[i].tempOccTime = 0.0;
                s[i].neighbors = NULL;
                s[i].nNeighbors = 0;
                s[i].rateSum = 0.0;
                s[i].updatedAt = 0.0;

            }

    return s;
}

/*
 * Distributes the electrons over random sites. Since all site energies
 * and positions were created randomly, simply the first nCarriers sites
 * can be set to occupied.m This can be modified to have a localized
 * source of electrons or anything else.
 */
Carrier *
MC_distributeCarriers (Site * sites)
{
    int i, k;
    Carrier *c;
    Carrier tmp;

    c = (Carrier *) malloc (sizeof (Carrier) * prms.ncarriers);

    for (k = 0; k < prms.ncarriers; ++k)
    {
        i = k;
        c[i].site = &sites[i];
        c[i].index = i;
        c[i].dx = 0.0;
        c[i].dy = 0.0;
        c[i].dz = 0.0;
        c[i].nFailedAttempts = 0;
        c[i].nHops = 0;
        c[i].occTime =
            (float) gsl_ran_exponential (prms.r, 1.0) / sites[i].rateSum;

        // now insert the carrier into the heap
        while (i > 0 && c[i].occTime < c[i / 2].occTime)
        {
            tmp = c[i];
            c[i] = c[i / 2];
            c[i / 2] = tmp;
            i /= 2;
        }

        c[i].site->carrier = &c[i];
        c[i].site->tempOccTime = 0.0001;
    }

    return c;
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
    SLE *head;

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
        temp->siteList = (SLE *) malloc (sizeof (SLE));
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
    SLE *sList, *tmp;
    Cell *cells = (Cell *) createCells (sites);

    for (l = 0; l < 100; ++l)
    {
        // this weirdness with 99 takes care of site numbers that
        // cannot be divided by 100.
        for (k = l * prms.nsites / 99; k < ((l + 1) * prms.nsites / 99); ++k)
            if (k < prms.nsites)
                setNeighbors (&sites[k], cells);

        if (serialOutput ())
        {
            printf ("\r\tInitializing...: \t\t%2d%%", (int) l);
            fflush (stdout);
        }
    }
    if (serialOutput ())
        printf (" Done.\n");

    // free cells
    for (l = 0; l <  prms.nx * prms.ny * prms.nz; ++l)
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
    int j, nSoftPairs = 0;
    bool ignore;
    SLE *ntemp = NULL, *neighbor = NULL;
    float rateb, rateab, ratebx = 0.0, rateSum;

    // find softpairs
    for (j = 0; j < prms.nsites; ++j)
    {
        neighbor = sites[j].neighbors;
        while (neighbor)
        {
            if (neighbor->rate / sites[j].rateSum > prms.softpairthreshold)
            {
                newSoftpair = (Softpair *) malloc (sizeof (Softpair));
                newSoftpair->i = &sites[j];
                newSoftpair->j = neighbor->s;
                newSoftpair->next = sp;
                sp = newSoftpair;
                nSoftPairs++;
            }
            neighbor = neighbor->next;
        }
    }

    if (nSoftPairs > 0)
    {
        // some output
        if (serialOutput ())
            printf ("\tRemoving Softpairs...:");

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
                neighbor = softpair->i->neighbors;
                while (neighbor)
                {
                    if (softpair->j->index == neighbor->s->index)
                    {
                        rateab = neighbor->rate;
                        rateb = neighbor->s->rateSum;
                        neighbor->rate = 0.0;
                        break;
                    }
                    neighbor = neighbor->next;
                }

                // update transition rates of all other neighbors
                rateSum = 0.0;
                neighbor = softpair->i->neighbors;
                while (neighbor)
                {
                    // find the partner
                    if (softpair->j->index != neighbor->s->index)
                    {
                        ratebx = 0.0;
                        ntemp = softpair->j->neighbors;
                        while (ntemp)
                        {
                            if (ntemp->s->index == neighbor->s->index)
                                ratebx = ntemp->rate;
                            ntemp = ntemp->next;
                        }
                        if (ratebx > 0)
                        {
                            neighbor->rate *= rateb;
                            neighbor->rate += ratebx * rateab;
                            neighbor->rate /= (rateb + rateab);
                        }
                        rateSum += neighbor->rate;
                    }
                    neighbor = neighbor->next;
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
        {
            neighbor = sites[j].neighbors;
            while (neighbor)
            {
                if (neighbor->rate / sites[j].rateSum > prms.softpairthreshold)
                    counter++;
                neighbor = neighbor->next;
            }
        }
        // some output
        if (serialOutput ())
            printf ("\t\tDone. %d Softpairs found.\n", counter);

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
    ssize_t i, k, l;
    Cell *c;
    SLE *siteList, *neighbors;
    Vector d;

    float rc2 = pow (prms.cutoff_radius, 2.0);

    neighbors = NULL;

    // make sure, cells are only used once
    int *used = malloc (27 * sizeof (int));
    int counter = 0, m, cont = 0;

    for (m = 0; m < 27; ++m)
        used[m] = -1;

    // now, find all the neighbors
    for (i = -1; i <= 1; i++)
        for (k = -1; k <= 1; k++)
            for (l = -1; l <= 1; l++)
            {

                c = getCell3D (cells, i + floor (s->x / prms.cutoff_radius),
                               floor (k + s->y / prms.cutoff_radius),
                               floor (l + s->z / prms.cutoff_radius));

                cont = 0;
                for (m = 0; m < 27; ++m)
                    if (used[m] == c->index)
                    {
                        cont = 1;
                        break;
                    }

                if (!cont)
                {
                    siteList = c->siteList;
                    used[counter++] = c->index;

                    while (siteList)
                    {
                        // if this is not the same Site as s and
                        // is within the cutoff radius, attach it 
                        // to the neighbor list and count one up.
                        d = distance (s, siteList->s);

                        if (s->index != siteList->s->index &&
                            pow (d.x, 2.0) + pow (d.y, 2.0) + pow (d.z,
                                                                   2.0) < rc2)
                        {

                            neighbors = (SLE *) malloc (sizeof (SLE));
                            neighbors->s = siteList->s;
                            neighbors->rate =
                                calcHoppingRate (*s, *siteList->s);
                            neighbors->dist = d;
                            neighbors->nTransitions = 0;
                            s->rateSum += neighbors->rate;

                            s->nNeighbors++;
                            neighbors->next = s->neighbors;
                            s->neighbors = neighbors;
                        }
                        siteList = siteList->next;
                    }
                }
            }

    free (used);

    // sort the neighbors according to the rate to save computation
    // time while simulating
    s->neighbors = sortNeighbors (s->neighbors);
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
    if (dE > 0)
        r *= exp (-1.0 * dE / prms.temperature);

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
    int lx, ly, lz;
    Vector vec;

    lx = prms.length_x / 2.;
    ly = prms.length_y / 2.;
    lz = prms.length_z / 2.;

    vec.x = j->x - i->x;
    vec.y = j->y - i->y;
    vec.z = j->z - i->z;

    if (vec.x > lx)
        vec.x -= prms.length_x;
    if (vec.x < -lx)
        vec.x += prms.length_x;

    if (vec.y > ly)
        vec.y -= prms.length_y;
    if (vec.y < -ly)
        vec.y += prms.length_y;

    if (vec.z > lz)
        vec.z -= prms.length_z;
    if (vec.z < -lz)
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

SLE *
sortNeighbors (SLE * list)
{
    if (!list || !list->next)
        return list;

    SLE *right = list, *temp = list, *last = list,
        *result = 0, *next = 0, *tail = 0;

    while (temp && temp->next)
    {
        last = right;
        right = right->next;
        temp = temp->next->next;
    }

    // Break the list in two.
    last->next = 0;

    // Recurse on the two smaller lists:
    list = sortNeighbors (list);
    right = sortNeighbors (right);

    // Merge:
    while (list || right)
    {
        // Take from empty lists, or compare:
        if (!right)
        {
            next = list;
            list = list->next;
        }
        else if (!list)
        {
            next = right;
            right = right->next;
        }
        else if (list->rate > right->rate)
        {
            next = list;
            list = list->next;
        }
        else
        {
            next = right;
            right = right->next;
        }

        if (!result)
            result = next;
        else
            tail->next = next;

        tail = next;
    }
    return result;
}
