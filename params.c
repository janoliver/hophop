#include "hop.h"

/*
 * Initialize all the params and make some rudimentary checks.
 */
void
generateParams (Params * prms, int argc, char **argv)
{

    // random number and timer stuff
    prms->loctime = localtime (&prms->curtime);
    gsl_rng_env_setup ();
    prms->T = gsl_rng_gfsr4;
    prms->r = gsl_rng_alloc (prms->T);

    // the gengetopt arguments
    struct gengetopt_args_info args;
    struct cmdline_parser_params *params;
    params = cmdline_parser_params_create ();
    prms->cmdlineargs = &args;

    // random seed
    prms->rseed = 0;
    if (args.rseed_given)
    {
        prms->rseed = args.rseed_arg;
    }

    // init command line parser and exit if anything went wrong
    if (cmdline_parser (argc, argv, &args) != 0)
    {
        printf ("Commandline error\n");
        exit (1);
    }

    // do we have to load a config file?
    params->initialize = 0;
    params->override = 0;
    if (args.conf_file_given &&
        cmdline_parser_config_file (args.conf_file_arg, &args, params) != 0)
    {
        printf ("Config file does not exist!\n");
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

    // check cutoff radius
    if (0 > args.rc_arg ||
        args.rc_arg > GSL_MIN (GSL_MIN (prms->length_y, prms->length_x),
                               prms->length_z))
    {
        printf ("Please choose a valid cut-off radius!\n");
        exit (1);
    }
    prms->cutoff_radius = args.rc_arg;

    // sigma
    if (0 >= args.sigma_arg)
    {
        printf ("Please choose a valid sigma!\n");
        exit (1);
    }
    prms->sigma = args.sigma_arg;

    // exponent
    if (0 >= args.exponent_arg)
    {
        printf ("Please choose a valid DOS exponent!\n");
        exit (1);
    }
    prms->exponent = args.exponent_arg;

    // loclength
    if (0 >= args.llength_arg || args.llength_arg > 2)
    {
        printf ("Please choose a valid localization length!\n");
        exit (1);
    }
    prms->loclength = args.llength_arg;

    // softpairthreshold
    if (0 >= args.softpairthreshold_arg)
    {
        printf ("Please choose a valid softpair threshold\n");
        exit (1);
    }
    prms->softpairthreshold = args.softpairthreshold_arg;

    // field
    prms->field = args.field_arg;

    // temperature
    if (0 >= args.temperature_arg)
    {
        printf ("Please choose a valid temperature!\n");
        exit (1);
    }
    prms->temperature = args.temperature_arg;

    // flags
    prms->gaussian = (args.gaussian_given) ? true : false;
    prms->accept_reject = (args.ar_given) ? true : false;
    prms->lattice = (args.lattice_given) ? true : false;
    prms->removesoftpairs = (args.removesoftpairs_given) ? true : false;
    prms->parallel = (args.parallel_given) ? true : false;
    prms->quiet = (args.quiet_given) ? true : false;
    prms->output_transitions = (args.transitions_given) ? true : false;
    prms->memreq = (args.memreq_given) ? true : false;

    // strings
    prms->output_folder = args.outputfolder_arg;
    prms->output_summary = args.summary_arg;
    prms->comment = args.comment_arg;

    //strcpy(prms->output_folder, args.outputfolder_arg);
    //strcpy(prms->output_summary, args.summary_arg);
    //strcpy(prms->comment, args.comment_arg);

    // simulation times
    if (0 >= args.relaxation_arg || 0 >= args.simulation_arg)
    {
        printf ("Please choose valid simulation times\n");
        exit (1);
    }
    prms->relaxation = args.relaxation_arg;
    prms->simulation = args.simulation_arg;

    // calculate number of cells
    prms->nx = ceil (prms->length_x / args.rc_arg);
    prms->ny = ceil (prms->length_y / args.rc_arg);
    prms->nz = ceil (prms->length_z / args.rc_arg);

    // check charge carriers
    if (args.ncarriers_arg >= prms->nsites)
    {
        printf ("Your system is overfilled with carriers.\n");
        exit (1);
    }
    prms->ncarriers = args.ncarriers_arg;

    // number of runs
    if (args.nruns_arg < 1)
        args.nruns_arg = 1;

    prms->number_runs = args.nruns_arg;

}

/*
 * Is a string argument given?
 */
bool
strArgGiven (char *arg)
{
    return (arg != NULL);
}

/*
 * Check if quiet is not set and we don't
 * use parallel computing atm.
 */
bool
serialOutput ()
{
    if (prms.quiet)
        return false;
    if (prms.parallel && prms.number_runs > 1)
        return false;
    return true;
}
