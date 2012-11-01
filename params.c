#include "hop.h"

/*
 * Initialize all the params and make some rudimentary checks.
 */
struct gengetopt_args_info args;

void
generateParams (Params * prms, int argc, char **argv)
{

    // the gengetopt arguments
    struct cmdline_parser_params *params;
    params = cmdline_parser_params_create ();
    prms->cmdlineargs = &args;

    // init command line parser and exit if anything went wrong
    if (cmdline_parser (argc, argv, &args) != 0)
    {
        output (O_FORCE, "Commandline error\n");
        exit (1);
    }

    // random number and timer stuff
    prms->loctime = localtime (&prms->curtime);
    gsl_rng_env_setup ();
    prms->T = gsl_rng_gfsr4;
    //prms->r = gsl_rng_alloc (prms->T);

    // random seed
    prms->rseed = 0;
    if (args.rseed_given)
    {
        prms->rseed = args.rseed_arg;
    }

    // do we have to load a config file?
    params->initialize = 0;
    params->override = 0;
    if (args.conf_file_given &&
        cmdline_parser_config_file (args.conf_file_arg, &args, params) != 0)
    {
        output (O_FORCE, "Config file does not exist!\n");
        exit (1);
    }

    // sample size stuff. The priority is as follows:
    // if x,y,z are given, use these.
    // if nsites is given, assume cubic sample and use it.
    // if length is given, assume cubic sample and use it.
    if (!(args.X_given && args.Y_given && args.Z_given))
    {
        if (args.nsites_given)
        {
            args.length_arg = floor (pow (args.nsites_arg, 1. / 3.));
            args.length_given = 1;
        }

        if (args.length_given)
        {
            args.X_arg = args.length_arg;
            args.Y_arg = args.length_arg;
            args.Z_arg = args.length_arg;
        }
    }

    prms->length_x = args.X_arg;
    prms->length_y = args.Y_arg;
    prms->length_z = args.Z_arg;

    prms->nsites = prms->length_x * prms->length_y * prms->length_z;

    // cutout stuff
    prms->cut_dos = false;
    if (args.cutoutenergy_given)
    {
        prms->cut_out_energy = args.cutoutenergy_arg;
        prms->cut_out_width = args.cutoutwidth_arg;
        prms->cut_dos = true;
    }

    // check cutoff radius
    if (0 > args.rc_arg ||
        args.rc_arg > GSL_MIN (GSL_MIN (prms->length_y, prms->length_x),
                               prms->length_z))
    {
        output (O_FORCE, "Please choose a valid cut-off radius!\n");
        exit (1);
    }
    prms->cutoff_radius = args.rc_arg;

    // exponent
    if (0 >= args.exponent_arg)
    {
        output (O_FORCE, "Please choose a valid DOS exponent!\n");
        exit (1);
    }
    prms->exponent = args.exponent_arg;

    // loclength
    if (0 >= args.llength_arg || args.llength_arg > 2)
    {
        output (O_FORCE, "Please choose a valid localization length!\n");
        exit (1);
    }
    prms->loclength = args.llength_arg;

    // softpairthreshold
    if (0 >= args.softpairthreshold_arg)
    {
        output (O_FORCE, "Please choose a valid softpair threshold\n");
        exit (1);
    }
    prms->softpairthreshold = args.softpairthreshold_arg;

    // field
    prms->field = args.field_arg;

    // temperature
    if (0 > args.temperature_arg)
    {
        output (O_FORCE, "Please choose a valid temperature!\n");
        exit (1);
    }
    prms->temperature = args.temperature_arg;

    // flags
    prms->gaussian = (args.gaussian_given) ? true : false;
    prms->lattice = (args.lattice_given) ? true : false;
    prms->removesoftpairs = (args.removesoftpairs_given) ? true : false;
    prms->parallel = (args.parallel_given) ? true : false;
    prms->quiet = (args.quiet_given) ? true : false;
    prms->output_transitions = (args.transitions_given) ? true : false;
    prms->memreq = (args.memreq_given) ? true : false;
    prms->balance_eq = (args.be_given) ? true : false;
    prms->meanfield = (args.meanfield_given) ? true : false;

    // balance equation parameters
    if (args.be_it_arg == 0)
        args.be_it_arg = prms->nsites - 1;
    prms->be_abs_tol = args.tol_abs_arg;
    prms->be_rel_tol = args.tol_rel_arg;
    prms->be_it = args.be_it_arg;
    prms->be_outer_it = args.be_oit_arg;

    // strings
    prms->output_folder = args.outputfolder_arg;
    prms->output_summary = args.summary_arg;
    prms->comment = args.comment_arg;

    // simulation times
    if (0 > args.relaxation_arg || 0 > args.simulation_arg)
    {
        output (O_FORCE, "Please choose valid simulation times\n");
        exit (1);
    }
    prms->relaxation = args.relaxation_arg;
    prms->simulation = args.simulation_arg;

    // calculate number of cells
    prms->nx = ceil (prms->length_x / prms->cutoff_radius);
    prms->ny = ceil (prms->length_y / prms->cutoff_radius);
    prms->nz = ceil (prms->length_z / prms->cutoff_radius);

    // check charge carriers
    if (args.ncarriers_arg >= prms->nsites)
    {
        output (O_FORCE, "Your system is overfilled with carriers.\n");
        exit (1);
    }
    prms->ncarriers = args.ncarriers_arg;

    // currently, the balance equations can only handle one single carrier.
    // therefore, if we use them, ncarriers is set to 1
    if (prms->balance_eq)
        prms->ncarriers = 1;

    // number of runs
    if (args.nruns_arg < 1)
        args.nruns_arg = 1;
    prms->number_runs = args.nruns_arg;

    // number of reruns (starting pos of the electron)
    if (args.nreruns_arg < 1)
        args.nreruns_arg = 1;
    prms->number_reruns = args.nreruns_arg;

    // threads
    if (args.nthreads_arg == 0 || args.nthreads_arg > omp_get_max_threads ())
        prms->nthreads = (GSL_MIN (omp_get_max_threads (), prms->number_runs));
    else
        prms->nthreads = (GSL_MIN (args.nthreads_arg, prms->number_runs));

    // if prms.nthreads == 1 for any reason, disable parallel computing
    if (prms->nthreads == 1)
        prms->parallel = false;

    // free memory
    free (params);

}

/*
 * Is a string argument given?
 */
bool
strArgGiven (char *arg)
{
    return (arg != NULL);
}
