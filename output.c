//  Created by Jan Oliver Oelerich

#include <sys/stat.h> 
#include "hop.h"
#include "mc.h"

void checkOutputFolder();

/**
    Writes the sites file with the following format:
    N    X   Y   Z
    x1   y1  z2  E1
    x2   y2  z2  E2
          .
          .
          .
    xn   yn  zn  En
    
    Filename: sites-config.dat
*/
void
writeSitesConfig(Site * sites)
{
    checkOutputFolder();
    
    FILE *file;
    int i;
    char fileName[128] = "";
    
    sprintf(fileName, "%s/sites-config.dat",args.outputfolder_arg);
    
    file = fopen(fileName,"w+");

    // write box information
    fprintf(file,"%8d\t%8d\t%8d\t%8d\n",args.nsites_arg,args.X_arg,args.Y_arg, args.Z_arg);

    // write site information
    for(i = 0;i<args.nsites_arg;++i)
    {
        fprintf(file,"%8.5f\t%8.5f\t%8.5f\t%8.5f\n",
            sites[i].x, sites[i].y, sites[i].z,sites[i].energy);
    }
    fclose(file);
    
    // some output
    if(!args.quiet_given)
        printf("Wrote site configuration to %s\n", fileName);
}

/**
    Writes the sites file with the following format:
    x1   y1  z2  E1  nvisited1   nvisitedUpwards1
    x2   y2  z2  E2  nvisited2   nvisitedUpwards2
          ...
    xn   yn  zn  En  nvisitedn   nvisitedUpwardsn
    
    Filename: sites.dat
*/
void
writeSites(Site * sites)
{
    checkOutputFolder();
    
    FILE *file;
    int i;
    char fileName[128] = "";
    
    sprintf(fileName, "%s/sites.dat",args.outputfolder_arg);
    
    file = fopen(fileName,"w+");

    // write site information
    for(i = 0;i<args.nsites_arg;++i)
    {
        fprintf(file,"%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8d\t%8d\n",
            sites[i].x, sites[i].y, sites[i].z,
            sites[i].energy, sites[i].visited,  sites[i].visitedUpward);
    }
    fclose(file);
    
    // some output
    if(!args.quiet_given)
        printf("Wrote site result information to %s\n", fileName);
}


/**
    Writes all the transitions to a datafile in the form

      index1 index2 E1 E2 NTransitions
 */
void
writeTransitions(Site * sites)
{
    checkOutputFolder();
    
    FILE *file;
    int i;
    char fileName[128] = "";
    
    sprintf(fileName, "%s/transitions.dat",args.outputfolder_arg);
    
    file = fopen(fileName,"w+");

    // write transitions information
    for(i = 0;i<args.nsites_arg;++i)
    {
        SLE * neighbor = sites[i].neighbors;
        while(neighbor)
        {
            if(neighbor->nTransitions > 0)
            {
                fprintf(file,"%d %d %8.5f %8.5f %8d\n",
                        sites[i].index, neighbor->s->index, sites[i].energy,
                        neighbor->s->energy, neighbor->nTransitions);
            }
            neighbor = neighbor->next;
        }
        
    }
    fclose(file);
    
    // some output
    if(!args.quiet_given)
        printf("Wrote transitions information to %s\n", fileName);
}


/**
    Writes the configuration file, which can be read
    by the program using the option --conf_file=/path/to/params.conf
    to re-use the settings.
    
    Filename: params.conf
*/
void
writeConfig()
{
    checkOutputFolder();
    
    char fileName[128] = "";
    
    sprintf(fileName, "%s/params.conf",args.outputfolder_arg);
    cmdline_parser_file_save(fileName, &args);
    
    // some output
    if(!args.quiet_given)
        printf("Wrote configuration file to %s\n", fileName);
}

/**
    Writes all the simulation results into one data file
    
    Filename: results.dat
*/
void 
writeResults(Results * res)
{
    checkOutputFolder();
    
    FILE *file, *file2;
    char fileName[128] = "";
    int buffer = 0;
    
    sprintf(fileName, "%s/results.dat",args.outputfolder_arg);
    
    // check for the header
    file2 = fopen(fileName,"r");
    if(file2 != NULL)
    {
        buffer = getc(file2);
        fclose(file2);
    }
    
    file = fopen(fileName,"a+");
    
    // write head
    if(buffer == 0)
    {
        fprintf(file, "#simul. time        ");
        fprintf(file, "mobility          ");
        fprintf(file, "diffusivity x/y   ");
        fprintf(file, "current density   ");
        fprintf(file, "fermi energy      ");
        fprintf(file, "random seed\n");
    }

    // write site information
    fprintf(file, "%f    ", res->simulationTime);
    fprintf(file, "%e      ", res->mobility);
    fprintf(file, "%e      ", res->diffusivity);
    fprintf(file, "%e      ", res->currentDensity);
    fprintf(file, "%f         ", res->fermiEnergy);
    fprintf(file, "%lu\n", RSEED);

    fclose(file);
    
    // some output
    if(!args.quiet_given)
        printf("\tWrote results to %s\n\n", fileName);
}

/**
    Checks if output folder argument is given. If not,
    create it's name and check if it exists. If not,
    create it.
*/
void
checkOutputFolder() 
{
    bool exists = false;
    
    // check if it exists
    struct stat st;
    if(stat(args.outputfolder_arg,&st) == 0)
        exists = true;
    
    if(!exists) 
    {
        char mkcmd[128];
        sprintf(mkcmd, "mkdir -p %s", args.outputfolder_arg);
        int ret = system(mkcmd);
        if(ret)
            printf("could not create output folder!\n");
    }
}

void
writeSummary(Results * res, Results * error)
{
    if(!args.summary_given)
        return;

    char fileName[1024];
    FILE *file, *file2;
    int buffer = 0;
    
    sprintf(fileName, "%s",args.summary_arg);
    
    // check for the header
    file2 = fopen(fileName,"r");
    if(file2 != NULL)
    {
        buffer = getc(file2);
        fclose(file2);
    }
    
    file = fopen(fileName,"a+");
    
    // write head
    if(buffer == 0)
    {
        fprintf(file, "#DOS exponent     ");
        fprintf(file, "System size       ");
        fprintf(file, "Number carriers   ");
        fprintf(file, "Localization len. ");
        fprintf(file, "Temperature       ");
        fprintf(file, "Electric field    ");
        fprintf(file, "Simul. time       ");
        fprintf(file, "Error simul. time ");
        fprintf(file, "Mobility          ");
        fprintf(file, "Error mobility    ");
        fprintf(file, "Diffusivity x/y   ");
        fprintf(file, "Error diffus. x/y ");
        fprintf(file, "Current density   ");
        fprintf(file, "Erro curr. dens.  ");
        fprintf(file, "Comment           ");
        fprintf(file, "\n");
    }

    // write site information
    fprintf(file, "%e      ", args.exponent_arg);
    fprintf(file, "%-10d        ", args.length_arg);
    fprintf(file, "%-10d        ", args.ncarriers_arg);
    fprintf(file, "%e      ", args.llength_arg);
    fprintf(file, "%e      ", args.temperature_arg);
    fprintf(file, "%e      ", args.field_arg);
    
    fprintf(file, "%e      ", res->simulationTime);
    fprintf(file, "%e      ", error->simulationTime);
    fprintf(file, "%e      ", res->mobility);
    fprintf(file, "%e      ", error->mobility);
    fprintf(file, "%e      ", res->diffusivity);
    fprintf(file, "%e      ", error->diffusivity);
    fprintf(file, "%e      ", res->currentDensity);
    fprintf(file, "%e      ", error->currentDensity);
    if(args.comment_given)
        fprintf(file, "%s", args.comment_arg);
    fprintf(file, "\n");
    
    fclose(file);
    
    // some output
    if(!args.quiet_given)
        printf("\nExtended summary file %s\n", fileName);
}
