//  Created by Jan Oliver Oelerich

#include "hop.h"

// the program parameters
Params prms;

void printResults (Results * results, Results * errors);
void printApplicationHeader ();
void printSettings ();
void checkApplicationSettings (int argc, char **argv);
void printEstimatedMemory ();

int
main (int argc, char **argv)
{
    // parse and check command line arguments
    // and generate the params struct which
    // is used all over the software
    generateParams (&prms, argc, argv);

    if (prms.memreq)
    {
        printEstimatedMemory ();
        return 0;
    }

    printApplicationHeader ();
    printSettings ();

    // start simulation
    int iRun;
    Results total[prms.number_runs], res, error;

    res.mobility = 0.0;
    res.diffusivity = 0.0;
    res.currentDensity = 0.0;
    res.simulationTime = 0.0;
    res.equilibrationEnergy = 0.0;

    if(prms.balance_eq) {
        beTest();
    }
    else
    {

        // set the number of threads
        omp_set_num_threads (GSL_MIN(omp_get_max_threads (),
                                     prms.number_runs));

#pragma omp parallel if(prms.parallel) shared(total) private(iRun)
        {
#pragma omp for schedule(dynamic)
            for (iRun = 1; iRun <= prms.number_runs; iRun++)
            {

                total[iRun - 1].mobility = 0.0;
                total[iRun - 1].diffusivity = 0.0;
                total[iRun - 1].currentDensity = 0.0;
                total[iRun - 1].simulationTime = 0.0;
                total[iRun - 1].equilibrationEnergy = 0.0;
                
                // here is where el magico happens
                MC_run (total, &iRun);
            }
        }
    
        // perform the averaging
        for (iRun = 0; iRun < prms.number_runs; iRun++)
        {
            res.mobility += total[iRun].mobility;
            res.diffusivity += total[iRun].diffusivity;
            res.currentDensity += total[iRun].currentDensity;
            res.simulationTime += total[iRun].simulationTime;
            res.equilibrationEnergy += total[iRun].equilibrationEnergy;
        }
        
    }

    
    res.mobility /= prms.number_runs;
    res.diffusivity /= prms.number_runs;
    res.currentDensity /= prms.number_runs;
    res.simulationTime /= prms.number_runs;
    res.equilibrationEnergy /= prms.number_runs;

    // find the errors
    for (iRun = 0; iRun < prms.number_runs; iRun++)
    {
        error.mobility = pow (res.mobility - total[iRun].mobility, 2);
        error.diffusivity = pow (res.diffusivity - total[iRun].diffusivity, 2);
        error.simulationTime = pow (res.simulationTime -
                                    total[iRun].simulationTime, 2);
        error.equilibrationEnergy = pow (res.equilibrationEnergy -
                                         total[iRun].equilibrationEnergy, 2);

    }
    error.mobility = sqrt (error.mobility / prms.number_runs);
    error.diffusivity = sqrt (error.diffusivity / prms.number_runs);
    error.currentDensity = sqrt (error.currentDensity / prms.number_runs);
    error.simulationTime = sqrt (error.simulationTime / prms.number_runs);
    error.equilibrationEnergy =
        sqrt (error.equilibrationEnergy / prms.number_runs);

    // summary
    if (strArgGiven (prms.output_summary))
        writeSummary (&res, &error);

    // output results to the command line
    printResults (&res, &error);

    gsl_rng_free (prms.r);

    return 0;
}

/*
 * This function prints out a useless header for the program.
 */
void
printApplicationHeader ()
{
    if (prms.quiet)
        return;
    printf ("\n\t###################################\n");
    printf ("\t#   Hopping Simulation            #\n");
    printf ("\t#   by Jan Oliver Oelerich, 2011  #\n");
    printf ("\t###################################\n");
}

/*
 * This function just prints out the parameters of the simulation.
 */
void
printSettings ()
{
    if (prms.quiet)
        return;
    // Settings output
    printf ("\nSettings:\n");
    printf ("\t3D sample size: \t\t%d x %d x %d\n",
            prms.length_x, prms.length_y, prms.length_z);
    printf ("\tDOS exponent: \t\t\tp = %1.1f\n", prms.exponent);
    printf ("\tDOS width: \t\t\ts = %2.4f\n", prms.sigma);
    printf ("\tLocalization length of sites: \ta = %2.4f\n", prms.loclength);
    printf ("\tTemperature: \t\t\tT = %2.4f\n", prms.temperature);
    printf ("\tField strength: \t\tF = %2.4f\n", prms.field);
    printf ("\tCut-off radius: \t\tr = %2.4f\n", prms.cutoff_radius);
    printf ("\tNumber of sites: \t\tN = %d\n", prms.nsites);
    printf ("\tNumber of carriers: \t\tn = %d\n", prms.ncarriers);
    printf ("\tHops of relaxation: \t\tR = %lu\n", prms.relaxation);
    printf ("\tHops of simulation: \t\tI = %lu\n", prms.simulation);
    if (!serialOutput () && !prms.quiet)
        printf ("\n");
}

/*
 * This function prints out the results of the simulation.
 */
void
printResults (Results * results, Results * error)
{
    if (prms.quiet)
        return;

    // Results output
    printf ("\nResults:\n");
    printf ("\tMobility in field-direction: \tu   = %e (+- %e)\n",
            results->mobility, error->mobility);
    printf ("\tDiffusivity perp. to field:  \tD   = %e (+- %e)\n",
            results->diffusivity, error->diffusivity);
    printf ("\tEinstein rel. perp. to field: \tD/u = %e\n",
            results->diffusivity / (results->mobility * prms.temperature));
    printf ("\tCurrent density (z-dir):  \tj   = %e (+- %e)\n",
            results->currentDensity, error->currentDensity);
    printf ("\tEquilibration Energy: \t\tE_i = %e\n",
            results->equilibrationEnergy);
    printf ("\tSimulated time: \t\tt   = %e\n\n", results->simulationTime);

}

void
printEstimatedMemory ()
{
    double mem = 0;

    // sites
    mem += prms.nsites * sizeof (Site);

    // carriers
    mem += prms.ncarriers * sizeof (Carrier);

    // neighbor lists
    mem +=
        pow (prms.cutoff_radius,
             3) * 4. / 3. * M_PI * prms.nsites * sizeof (SLE);

    // parallelization
    if (prms.parallel && prms.number_runs >= omp_get_max_threads ())
        mem *= omp_get_max_threads ();
    else if (prms.parallel)
        mem *= prms.number_runs;

    // results
    mem += (prms.number_runs + 2) * sizeof (Results);

    // some offset that shouldn't scale with the system...
    mem += 100 * 1024 * 1024;

    printf ("Estimated memory usage: %5.2f MB\n", mem / (1024 * 1024));
    return;
}
