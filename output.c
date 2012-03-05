//  Created by Jan Oliver Oelerich

#include <sys/stat.h>
#include "hop.h"

void checkOutputFolder (int iRun);

/*
 * Writes the sites file with the following format:
 * N    X   Y   Z
 * x1   y1  z2  E1
 * x2   y2  z2  E2
 * .
 * .
 * .
 * xn   yn  zn  En
 * 
 * Filename: sites-config.dat
 */
void
writeSitesConfig (Site * sites, int iRun)
{
    checkOutputFolder (iRun);

    FILE *file;
    int i;
    char fileName[128] = "";

    sprintf (fileName, "%s/%d/sites-config.dat", prms.output_folder, iRun);

    file = fopen (fileName, "w+");

    // write box information
    fprintf (file, "%8d\t%8d\t%8d\t%8d\n", prms.nsites, prms.length_x,
             prms.length_y, prms.length_z);

    // write site information
    for (i = 0; i < prms.nsites; ++i)
    {
        fprintf (file, "%8.5f\t%8.5f\t%8.5f\t%8.5f\n",
                 sites[i].x, sites[i].y, sites[i].z, sites[i].energy);
    }
    fclose (file);

    // some output
    if (serialOutput ())
        printf ("\tWrote site configuration to \t\t%s\n", fileName);
}

/*
 * Writes the sites file with the following format:
 * x1   y1  z2  E1  nvisited1   nvisitedUpwards1
 * x2   y2  z2  E2  nvisited2   nvisitedUpwards2
 * ...
 * xn   yn  zn  En  nvisitedn   nvisitedUpwardsn
 * 
 * Filename: sites.dat
 */
void
writeSites (Site * sites, int iRun)
{
    checkOutputFolder (iRun);

    FILE *file;
    int i;
    char fileName[128] = "";

    sprintf (fileName, "%s/%d/sites.dat", prms.output_folder, iRun);

    file = fopen (fileName, "w+");

    // write site information
    for (i = 0; i < prms.nsites; ++i)
    {
        fprintf (file, "%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8d\t%8d\n",
                 sites[i].x, sites[i].y, sites[i].z,
                 sites[i].energy, sites[i].visited, sites[i].visitedUpward);
    }
    fclose (file);

    // some output
    if (serialOutput ())
        printf ("\tWrote site result information to \t%s\n", fileName);
}


/*
 * Writes all the transitions to a datafile in the form
 * index1 index2 E1 E2 NTransitions
 */
void
writeTransitions (Site * sites, int iRun)
{
    checkOutputFolder (iRun);

    FILE *file;
    int i;
    char fileName[128] = "";

    sprintf (fileName, "%s/%d/transitions.dat", prms.output_folder, iRun);

    file = fopen (fileName, "w+");

    // write transitions information
    for (i = 0; i < prms.nsites; ++i)
    {
        SLE *neighbor = sites[i].neighbors;
        while (neighbor)
        {
            if (neighbor->nTransitions > 0)
            {
                fprintf (file, "%d %d %8.5f %8.5f %8d\n",
                         sites[i].index, neighbor->s->index, sites[i].energy,
                         neighbor->s->energy, neighbor->nTransitions);
            }
            neighbor = neighbor->next;
        }

    }
    fclose (file);

    // some output
    if (serialOutput ())
        printf ("\tWrote transitions information to \t%s\n", fileName);
}


/*
 * Writes the configuration file, which can be read by the program using
 * the option --conf_file=/path/to/params.conf to re-use the settings.
 * 
 * Filename: params.conf
 */
void
writeConfig (int iRun)
{
    checkOutputFolder (iRun);

    char fileName[128] = "";

    sprintf (fileName, "%s/params.conf", prms.output_folder);
    cmdline_parser_file_save (fileName, prms.cmdlineargs);

    // some output
    if (serialOutput ())
        printf ("\tWrote configuration file to \t\t%s\n", fileName);
}

/*
 * Writes all the simulation results into one data file
 * 
 * Filename: results.dat
 */
void
writeResults (Results * res, int iRun)
{
    checkOutputFolder (iRun);

    FILE *file, *file2;
    char fileName[128] = "";
    int buffer = 0;

    sprintf (fileName, "%s/results.dat", prms.output_folder);

    // check for the header
    file2 = fopen (fileName, "r");
    if (file2 != NULL)
    {
        buffer = getc (file2);
        fclose (file2);
    }

    file = fopen (fileName, "a+");

    // write head
    if (buffer == 0)
    {
        fprintf (file, "#simul. time      ");
        fprintf (file, "mobility          ");
        fprintf (file, "diffusivity x/y   ");
        fprintf (file, "current density   ");
        fprintf (file, "Equilibration en. ");
        fprintf (file, "random seed\n");
    }

    // write site information
    if(prms.balance_eq)
    {
        fprintf (file, "%-18e", res->simulationTime);
        fprintf (file, "%-+18e", res->mobility);
        fprintf (file, "%-+18e", res->diffusivity);
        fprintf (file, "%-+18e", res->currentDensity);
        fprintf (file, "%-+18e", res->equilibrationEnergy);
        fprintf (file, "%lu\n", prms.rseed_used);
    }
    else
    {
        fprintf (file, "%-+18e", res->simulationTime);
        fprintf (file, "%-+18e", res->mobility);
        fprintf (file, "%-+18e", res->diffusivity);
        fprintf (file, "%-+18e", res->currentDensity);
        fprintf (file, "%-+18e", res->equilibrationEnergy);
        fprintf (file, "%lu\n", prms.rseed_used);
    }
    
    fclose (file);

    // some output
    if (serialOutput ())
        printf ("\n\tWrote results to \t\t\t%s\n", fileName);
}

/*
 * Checks if output folder argument is given. If not, create it's name and
 * check if it exists. If not, create it.
 */
void
checkOutputFolder (int iRun)
{
    // check if realization folder exists
    char realfolder[200];
    sprintf (realfolder, "mkdir -p %s/%d", prms.output_folder, iRun);
    int ret = system (realfolder);
    if (ret)
        printf ("could not create output realization folder!\n");

}

void
writeSummary (Results * res, Results * error)
{
    if (!strArgGiven (prms.output_summary))
        return;

    char fileName[1024];
    FILE *file, *file2;
    int buffer = 0;

    sprintf (fileName, "%s", prms.output_summary);

    // check for the header
    file2 = fopen (fileName, "r");
    if (file2 != NULL)
    {
        buffer = getc (file2);
        fclose (file2);
    }

    file = fopen (fileName, "a+");

    // write head
    if (buffer == 0)
    {
        fprintf (file, "#DOS exponent     ");
        fprintf (file, "System size       ");
        fprintf (file, "Number carriers   ");
        fprintf (file, "Localization len. ");
        fprintf (file, "Temperature       ");
        fprintf (file, "Electric field    ");
        fprintf (file, "Simul. time       ");
        fprintf (file, "Mobility          ");
        fprintf (file, "Error mobility    ");
        fprintf (file, "Diffusivity x/y   ");
        fprintf (file, "Error diffus. x/y ");
        fprintf (file, "Current density   ");
        fprintf (file, "Erro curr. dens.  ");
        fprintf (file, "Equilibration en. ");
        fprintf (file, "Error Eq. en.     ");
        fprintf (file, "Average energy    ");
        fprintf (file, "Error avg. en.    ");
        fprintf (file, "Mode              ");
        fprintf (file, "Comment           ");
        fprintf (file, "\n");
    }

    // write site information
    fprintf (file, "%-+18e", prms.exponent);
    fprintf (file, "%-18d", prms.length_x);
    fprintf (file, "%-18d", prms.ncarriers);
    fprintf (file, "%-+18e", prms.loclength);
    fprintf (file, "%-+18e", prms.temperature);
    fprintf (file, "%-+18e", prms.field);

    fprintf (file, "%-+18e", res->simulationTime);
    fprintf (file, "%-+18e", res->mobility);
    fprintf (file, "%-+18e", error->mobility);
    fprintf (file, "%-+18e", res->diffusivity);
    fprintf (file, "%-+18e", error->diffusivity);
    fprintf (file, "%-+18e", res->currentDensity);
    fprintf (file, "%-+18e", error->currentDensity);
    fprintf (file, "%-+18e", res->equilibrationEnergy);
    fprintf (file, "%-+18e", error->equilibrationEnergy);
    fprintf (file, "%-+18e", res->avgenergy);
    fprintf (file, "%-+18e", error->avgenergy);
    fprintf (file, "%-18s",
             prms.balance_eq ? "Balance eq." : "MC Simulation");
    if (strArgGiven (prms.comment))
        fprintf (file, "%s", prms.comment);
    fprintf (file, "\n");

    fclose (file);

    // some output
    if (serialOutput ())
        printf ("\nExtended summary file %s\n", fileName);
}
