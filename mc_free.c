//  Created by Jan Oliver Oelerich

#include "mc.h"

void
MC_freeSites(Site * sites)
{
    SLE * neighbor, * tmp;
    int i;
    for(i = 0; i < args.nsites_arg; ++i)
    {
        // free neighbor memory
        neighbor = sites[i].neighbors;
        while(neighbor)
        {
            tmp = neighbor->next;
            free(neighbor);
            neighbor = tmp;
        }
    }
    free((void *)sites);

}

void
MC_freeCarriers(Carrier * carriers)
{
    free((void *)carriers);
}
