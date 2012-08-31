//  Created by Jan Oliver Oelerich

#include <sys/stat.h>
#include "hop.h"

void checkOutputFolder (RunParams * runprms);


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
writeSitesConfig (Site * sites, RunParams * runprms)
{
    checkOutputFolder (runprms);

    FILE *file;
    int i;
    char fileName[128] = "";

    sprintf (fileName, "%s/%d/sites-config.dat", prms.output_folder,
             runprms->iRun);

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
    output (O_SERIAL, "\tWrote site configuration to \t\t%s\n", fileName);
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
writeSites (Site * sites, RunParams * runprms)
{
    checkOutputFolder (runprms);

    FILE *file;
    int i;
    char fileName[128] = "";

    sprintf (fileName, "%s/%d/sites.dat", prms.output_folder, runprms->iRun);

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
    output (O_SERIAL, "\tWrote site result information to \t%s\n", fileName);
}


/*
 * Writes all the transitions to a datafile in the form
 * index1 index2 E1 E2 NTransitions
 */
void
writeTransitions (Site * sites, RunParams * runprms)
{
    checkOutputFolder (runprms);

    FILE *file;
    int i, j;
    char fileName[128] = "";
    SLE *neighbor;

    sprintf (fileName, "%s/%d/transitions.dat", prms.output_folder,
             runprms->iRun);

    file = fopen (fileName, "w+");

    // write transitions information
    for (i = 0; i < prms.nsites; ++i)
    {
        for (j = 0; j < sites[i].nNeighbors; ++j)
        {
            neighbor = &(sites[i].neighbors[j]);
            if (neighbor->nTransitions > 0)
            {
                fprintf (file, "%d %d %8.5f %8.5f %8d\n",
                         sites[i].index, neighbor->s->index, sites[i].energy,
                         neighbor->s->energy, neighbor->nTransitions);
            }
        }

    }
    fclose (file);

    // some output
    output (O_SERIAL, "\tWrote transitions information to \t%s\n", fileName);
}


/*
 * Writes the configuration file, which can be read by the program using
 * the option --conf_file=/path/to/params.conf to re-use the settings.
 * 
 * Filename: params.conf
 */
void
writeConfig (RunParams * runprms)
{
    checkOutputFolder (runprms);

    char fileName[128] = "";

    sprintf (fileName, "%s/params.conf", prms.output_folder);
    cmdline_parser_file_save (fileName, prms.cmdlineargs);

    // some output
    output (O_SERIAL, "\tWrote configuration file to \t\t%s\n", fileName);
}

/*
 * Writes all the simulation results into one data file
 * 
 * Filename: results.dat
 */
void
writeResults (Results * res, RunParams * runprms)
{
    checkOutputFolder (runprms);

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
        if (prms.balance_eq)
        {
            fprintf (file, "#mobility         ");
            fprintf (file, "random seed\n");
        }
        else
        {
            fprintf (file, "#simul. time      ");
            fprintf (file, "mobility          ");
            fprintf (file, "diffusivity x/y   ");
            fprintf (file, "current density   ");
            fprintf (file, "Equilibration en. ");
            fprintf (file, "random seed\n");
        }
    }

    // write site information
    if (prms.balance_eq)
    {
        fprintf (file, "%-+18e", res->mobility.values[runprms->iRun]);
        fprintf (file, "%lu\n", runprms->rseed_used);
    }
    else
    {
        fprintf (file, "%-18e", runprms->simulationTime);
        fprintf (file, "%-+18e", res->mobility.values[runprms->iRun]);
        fprintf (file, "%-+18e", res->diffusivity.values[runprms->iRun]);
        fprintf (file, "%-+18e", res->currentDensity.values[runprms->iRun]);
        fprintf (file, "%-+18e",
                 res->equilibrationEnergy.values[runprms->iRun]);
        fprintf (file, "%lu\n", runprms->rseed_used);
    }

    fclose (file);

    // some output
    output (O_SERIAL, "\n\tWrote results to \t\t\t%s\n", fileName);
}

/*
 * Checks if output folder argument is given. If not, create it's name and
 * check if it exists. If not, create it.
 */
void
checkOutputFolder (RunParams * runprms)
{
    // check if realization folder exists
    char realfolder[200];
    sprintf (realfolder, "mkdir -p %s/%d", prms.output_folder, runprms->iRun);
    int ret = system (realfolder);
    if (ret)
        output (O_FORCE, "could not create output realization folder!\n");

}

void
writeSummary (Results * res)
{
    if (!strArgGiven (prms.output_summary))
        return;

    char fileName[1024];
    FILE *file, *file2;
    int buffer = 0;

    // build the current time string
    char timestring[20]; 
    time_t now;
    time(&now);
    struct tm *ti;
    ti = localtime(&now);
    sprintf(timestring, "%04d-%02d-%02d_%02d:%02d:%02d",
            ti->tm_year+1900, ti->tm_mon, ti->tm_mday, ti->tm_hour,
            ti->tm_min, ti->tm_sec);
 
    sprintf (fileName, "%s", prms.output_summary);

    // check for the header
    file2 = fopen (fileName, "r");
    if (file2 != NULL)
    {
        buffer = getc (file2);
        fclose (file2);
    }

    file = fopen (fileName, "a+");

    // write header
    if (buffer == 0)
    {
        fprintf (file, "%s%s %s\n", "#", PKG_NAME, PKG_VERSION);
        
        fprintf (file, "%-20s", "#mode");
        fprintf (file, "%-20s", "dos_exponent");
        fprintf (file, "%-20s", "system_length");
        fprintf (file, "%-20s", "number_sites");
        fprintf (file, "%-20s", "number_carriers");
        fprintf (file, "%-20s", "localization_length");
        fprintf (file, "%-20s", "temperature");
        fprintf (file, "%-20s", "electric_field");
        fprintf (file, "%-20s", "number_runs");
        fprintf (file, "%-20s", "number_reruns");
        fprintf (file, "%-20s", "number_relaxation");
        fprintf (file, "%-20s", "number_simulation");
        fprintf (file, "%-20s", "cutout_energy");
        fprintf (file, "%-20s", "cutout_energy_err");
        fprintf (file, "%-20s", "simulation_time");
        fprintf (file, "%-20s", "simulation_time_err");
        fprintf (file, "%-20s", "mobility");
        fprintf (file, "%-20s", "mobility_err");
        fprintf (file, "%-20s", "diffusivity");
        fprintf (file, "%-20s", "diffusivity_err");
        fprintf (file, "%-20s", "einstein_rel");
        fprintf (file, "%-20s", "einstein_rel_err");
        fprintf (file, "%-20s", "current_dens");
        fprintf (file, "%-20s", "current_dens_err");
        fprintf (file, "%-20s", "eq_energy");
        fprintf (file, "%-20s", "eq_energy_err");
        fprintf (file, "%-20s", "avg_energy");
        fprintf (file, "%-20s", "avg_energy_err");
        fprintf (file, "%-20s", "random_seed");
        fprintf (file, "%-20s", "finish_time");
        fprintf (file, "%-20s", "comment");
        fprintf (file, "\n");

        fprintf (file, "%-20s", "#str");
        fprintf (file, "%-20s", "int");
        fprintf (file, "%-20s", "int");
        fprintf (file, "%-20s", "long");
        fprintf (file, "%-20s", "long");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "int");
        fprintf (file, "%-20s", "int");
        fprintf (file, "%-20s", "long");
        fprintf (file, "%-20s", "long");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "float");
        fprintf (file, "%-20s", "long");
        fprintf (file, "%-20s", "str");
        fprintf (file, "%-20s", "str");
        fprintf (file, "\n");
    }

    // write site information
    fprintf (file, "%-20s", prms.balance_eq ? "be" : "mc");
    fprintf (file, "%-+20e", prms.exponent);
    fprintf (file, "%-20d", prms.length_x);
    fprintf (file, "%-20d", prms.nsites);
    fprintf (file, "%-20d", prms.ncarriers);
    fprintf (file, "%-+20e", prms.loclength);
    fprintf (file, "%-+20e", prms.temperature);
    fprintf (file, "%-+20e", prms.field);
    fprintf (file, "%-20d", prms.number_runs);
    fprintf (file, "%-20d", prms.number_reruns);
    fprintf (file, "%-20lu", prms.relaxation);
    fprintf (file, "%-20lu", prms.simulation);
    fprintf (file, "%-+20e", prms.cut_dos ? prms.cut_out_energy : 0);
    fprintf (file, "%-+20e", prms.cut_dos ? prms.cut_out_width : 0);
    fprintf (file, "%-+20e", res->simulationTime.avg);
    fprintf (file, "%-+20e", res->simulationTime.err);
    fprintf (file, "%-+20e", res->mobility.avg);
    fprintf (file, "%-+20e", res->mobility.err);
    fprintf (file, "%-+20e", res->diffusivity.avg);
    fprintf (file, "%-+20e", res->diffusivity.err);
    fprintf (file, "%-+20e", res->einsteinrelation.avg);
    fprintf (file, "%-+20e", res->einsteinrelation.err);
    fprintf (file, "%-+20e", res->currentDensity.avg);
    fprintf (file, "%-+20e", res->currentDensity.err);
    fprintf (file, "%-+20e", res->equilibrationEnergy.avg);
    fprintf (file, "%-+20e", res->equilibrationEnergy.err);
    fprintf (file, "%-+20e", res->avgenergy.avg);
    fprintf (file, "%-+20e", res->avgenergy.err);
    fprintf (file, "%-20lu", prms.rseed);
    fprintf (file, "%-20s", timestring);
    
    if (strArgGiven (prms.comment))
        fprintf (file, "%s", prms.comment);
    fprintf (file, "\n");

    fclose (file);

    // some output
    output (O_SERIAL, "\nExtended summary file %s\n", fileName);
}
