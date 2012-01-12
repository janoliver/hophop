//  Created by Jan Oliver Oelerich

#include "hop.h"

struct gengetopt_args_info args;
unsigned long RSEED;
time_t curtime;
const gsl_rng_type * T;
gsl_rng * r;
int nx, ny, nz;

void printResults(Results * results, Results * errors);
void printApplicationHeader();
void printSettings();
void checkApplicationSettings(int argc, char **argv);
void printEstimatedMemory();

int
main(int argc, char **argv)
{
    loctime = localtime(&curtime);

    checkApplicationSettings(argc, argv);

    if(args.memreq_given)
    {
        printEstimatedMemory();
        return 0;
    }
    
    printApplicationHeader();
    printSettings();

    // start simulation
    int iRun;
    Results total[args.nruns_arg], res, error;

    res.mobility         = 0.0;
    res.diffusivity      = 0.0;
    res.fermiEnergy      = 0.0;
    res.transportEnergy  = 0.0;
    res.currentDensity   = 0.0;
    res.simulationTime   = 0.0;
    res.equilibrationEnergy = 0.0;    
    
    // random initialization
    gsl_rng_env_setup();
    T = gsl_rng_gfsr4;
    r = gsl_rng_alloc(T);
    
    // averaging loop
    int nThreads = 0;
    if(args.nruns_given && args.nruns_arg <= omp_get_max_threads())
        nThreads = args.nruns_arg;
    else
        nThreads = omp_get_max_threads();
    
    omp_set_num_threads(nThreads);
#pragma omp parallel if(args.parallel_given) shared(total) private(iRun)
    {
#pragma omp for schedule(dynamic)
        for(iRun = 1; iRun <= args.nruns_arg; iRun++)
        {
        
            total[iRun-1].mobility         = 0.0;
            total[iRun-1].diffusivity      = 0.0;
            total[iRun-1].fermiEnergy      = 0.0;
            total[iRun-1].transportEnergy  = 0.0;
            total[iRun-1].currentDensity   = 0.0;
            total[iRun-1].simulationTime   = 0.0;
            total[iRun-1].equilibrationEnergy = 0.0;
            
            // here is where el magico happens
            MC_run(total, &iRun);
        }
    }
    // perform the averaging
    for(iRun = 0; iRun < args.nruns_arg; iRun++)
    {
        res.mobility        += total[iRun].mobility;
        res.diffusivity     += total[iRun].diffusivity;
        res.fermiEnergy     += total[iRun].fermiEnergy;
        res.transportEnergy += total[iRun].transportEnergy;
        res.currentDensity  += total[iRun].currentDensity;
        res.simulationTime  += total[iRun].simulationTime;
        res.equilibrationEnergy += total[iRun].equilibrationEnergy;
    }
    res.mobility        /= args.nruns_arg;
    res.diffusivity     /= args.nruns_arg;
    res.fermiEnergy     /= args.nruns_arg;
    res.transportEnergy /= args.nruns_arg;
    res.currentDensity  /= args.nruns_arg;
    res.simulationTime  /= args.nruns_arg;
    res.equilibrationEnergy /= args.nruns_arg;

    // find the errors
    for(iRun = 0; iRun < args.nruns_arg; iRun++)
    {
        error.mobility        = pow(res.mobility - total[iRun].mobility, 2);
        error.diffusivity     = pow(res.diffusivity -
                                    total[iRun].diffusivity, 2);
        error.fermiEnergy     = pow(res.fermiEnergy -
                                    total[iRun].fermiEnergy, 2);
        error.transportEnergy = pow(res.transportEnergy -
                                    total[iRun].transportEnergy, 2);
        error.simulationTime  = pow(res.simulationTime -
                                    total[iRun].simulationTime, 2);
        error.equilibrationEnergy  = pow(res.equilibrationEnergy -
                                         total[iRun].equilibrationEnergy, 2);
        
    }
    error.mobility        = sqrt(error.mobility/args.nruns_arg);
    error.diffusivity     = sqrt(error.diffusivity/args.nruns_arg);
    error.fermiEnergy     = sqrt(error.fermiEnergy/args.nruns_arg);
    error.transportEnergy = sqrt(error.transportEnergy/args.nruns_arg);
    error.currentDensity  = sqrt(error.currentDensity/args.nruns_arg);
    error.simulationTime  = sqrt(error.simulationTime/args.nruns_arg);
    error.equilibrationEnergy = sqrt(error.equilibrationEnergy/args.nruns_arg);
    
    // summary
    if(args.summary_given)
        writeSummary(&res, &error);

    // output results to the command line
    printResults(&res, &error);

    gsl_rng_free(r);

    return 0;
}

/**
    This function checks and processes application settings.
*/
void
checkApplicationSettings(int argc, char **argv)
{
    // initialize the parameters structure
    struct cmdline_parser_params *params;
    params = cmdline_parser_params_create();

    // init command line parser and exit if anything went wrong
    if(cmdline_parser(argc, argv, &args) != 0)
    {
        printf("Commandline error\n");
        exit(1);
    }

    // config file
    params->initialize = 0;
    params->override   = 0;
    if(args.conf_file_given &&
       cmdline_parser_config_file(args.conf_file_arg, &args, params) != 0)
    {
        printf("Config file does not exist!\n");
        exit(1);
    }

    // set lenths if the length option is given
    if(args.length_given)
    {
        args.X_arg = args.length_arg;
        args.Y_arg = args.length_arg;
        args.Z_arg = args.length_arg;
    }

    // set sample size according to site concentration
    // and the other way around
    if((args.length_given || (args.X_given && args.Y_given &&
                              args.Z_given))
       && !args.nsites_given)
    {
        args.nsites_arg = args.X_arg * args.Y_arg * args.Z_arg;
    }
    if(!(args.length_given || (args.X_given && args.Y_given && args.Z_given))
       && args.nsites_given)
    {
        int sitesPerDirection = pow(args.nsites_arg,1./3.);
        args.X_arg = sitesPerDirection;
        args.Y_arg = sitesPerDirection;
        args.Z_arg = sitesPerDirection;
    }

    // calculate number of cells
    nx = ceil(args.X_arg / args.rc_arg);
    ny = ceil(args.Y_arg / args.rc_arg);
    nz = ceil(args.Z_arg / args.rc_arg);

    // error checking..
    if(args.ncarriers_arg >= args.nsites_arg)
    {
        printf("Your system is overfilled with carriers.\n");
        exit(1);
    }

    if(args.nruns_arg < 1)
        args.nruns_arg = 1;
}

/**
   This function prints out a useless header for the program.
*/
void
printApplicationHeader()
{
    if(args.quiet_given)
        return;
    printf("\n\t###################################\n");
    printf("\t#   Hopping Simulation            #\n");
    printf("\t#   by Jan Oliver Oelerich, 2011  #\n");
    printf("\t###################################\n");
}

/**
   This function just prints out the parameters of the simulation.
*/
void
printSettings()
{
    if(args.quiet_given)
        return;
    // Settings output
    printf("\nSettings:\n");
    printf("\t3D sample size: \t\t%d x %d x %d\n",
           args.X_arg, args.Y_arg, args.Z_arg);
    printf("\tDOS exponent: \t\t\tp = %1.1f\n",
           args.exponent_arg);
    printf("\tDOS width: \t\t\ts = %2.4f\n",
           args.sigma_arg);
    printf("\tLocalization length of sites: \ta = %2.4f\n",
           args.llength_arg);
    printf("\tTemperature: \t\t\tT = %2.4f\n",
           args.temperature_arg);
    printf("\tField strength: \t\tF = %2.4f\n",
           args.field_arg);
    printf("\tCut-off radius: \t\tr = %2.4f\n",
           args.rc_arg);
    printf("\tNumber of sites: \t\tN = %d\n",
           args.nsites_arg);
    printf("\tNumber of carriers: \t\tn = %d\n",
           args.ncarriers_arg);
    printf("\tHops of relaxation: \t\tR = %lu\n",
           args.relaxation_arg);
    printf("\tHops of simulation: \t\tI = %lu\n",
           args.simulation_arg);
    if(!args.quiet_given && args.parallel_given && args.nruns_arg > 1)
        printf("\n");
}

/**
   This function prints out the results of the simulation.
*/
void
printResults(Results * results, Results * error)
{
    if(args.quiet_given)
        return;
    
    // Results output
    printf("\nResults:\n");
    printf("\tMobility in field-direction: \tu   = %e (+- %e)\n",
           results->mobility, error->mobility);
    printf("\tDiffusivity perp. to field:  \tD   = %e (+- %e)\n",
           results->diffusivity, error->diffusivity);
    printf("\tEinstein rel. perp. to field: \tD/u = %e\n",
           results->diffusivity / (results->mobility * args.temperature_arg));
    printf("\tCurrent density (z-dir):  \tj   = %e (+- %e)\n",
           results->currentDensity, error->currentDensity);
    printf("\tFermi energy:  \t\t\tE_f = %e\n",
      results->fermiEnergy);
    printf("\tEquilibration Energy: \t\tE_i = %e\n",
           results->equilibrationEnergy);
    printf("\tSimulated time: \t\tt   = %e\n\n",
           results->simulationTime);

}

void
printEstimatedMemory()
{
    double mem = 0;
    
    // sites
    mem += args.nsites_arg * sizeof(Site);

    // carriers
    mem += args.ncarriers_arg * sizeof(Carrier);

    // neighbor lists
    mem += pow(args.rc_arg, 3) * 4. / 3. * M_PI * args.nsites_arg * sizeof(SLE);

    // parallelization
    if(args.parallel_given && args.nruns_given >= omp_get_max_threads())
        mem *= omp_get_max_threads();
    else if(args.parallel_given)
        mem *= args.nruns_arg;

    // results
    mem += (args.nruns_arg+2) * sizeof(Results);

    // some offset that shouldn't scale with the system...
    mem += 100*1024*1024;

    printf("Estimated memory usage: %5.2f MB\n",mem/(1024*1024));
    return;
}
