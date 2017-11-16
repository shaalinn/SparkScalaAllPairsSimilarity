/*!
 \file  cmdline.c
 \brief This file contains functions for parsing command-line arguments

 \author David C. Anastasiu
 */
#include "includes.h"


/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
	{"t",                 1,      0,      CMD_SIMTHRESHOLD},
	{"eps",               1,      0,      CMD_EPSILON},
    {"sim",               1,      0,      CMD_SIM},
	{"nsk",               1,      0,      CMD_NSKETCHES},
	{"fldelta",           1,      0,      CMD_FLDELTA},
	{"nim",               0,      0,      CMD_NIM},
	{"verb",              1,      0,      CMD_VERBOSITY},
	{"v",                 0,      0,      CMD_VERSION},
	{"pr",                0,      0,      CMD_PR},
	{"prmin",             1,      0,      CMD_PRMINLEN},
	{"prmax",             1,      0,      CMD_PRMAXLEN},
	{"pc",                0,      0,      CMD_PC},
	{"pcmin",             1,      0,      CMD_PCMINLEN},
	{"pcmax",             1,      0,      CMD_PCMAXLEN},
	{"fmtRead",           1,      0,      CMD_FMT_READ},
	{"readZidx",          0,      0,      CMD_FMT_READ_NUM},
	{"readVals",          1,      0,      CMD_READ_VALS},
	{"fmtWrite",          1,      0,      CMD_FMT_WRITE},
	{"writeZidx",         0,      0,      CMD_FMT_WRITE_NUM},
	{"writeVals",         1,      0,      CMD_WRITE_VALS},
	{"compactCols",       0,      0,      CMD_COMPACT_COLS},
	{"compactRows",       0,      0,      CMD_COMPACT_ROWS},
	{"seed",              1,      0,      CMD_SEED},
	{"scale",             0,      0,      CMD_SCALE},
	{"norm",              1,      0,      CMD_NORM},

	{"help",              0,      0,      CMD_HELP},
	{"h",                 0,      0,      CMD_HELP},
	{0,                   0,      0,      0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: " PROGRAM_NAME " [options] mode input-file [output-file]",
" ",
" mode:",
"  ij      IdxJoin - Full sparse dot-product with lesser id docs.",
"  ap      AllPairs (Cosine only)",
"  mmj     MMJoin",
"  mkj     MK-Join (Tanimoto only)",
"  mkj2    MK-Join with tighter bounds from ACIIDS 2014 paper",
"  l2ap    L2-Norm AllPairs (L2AP)",
"  tapnn   Tanimoto extension of L2AP",
"  l2ap-c, tapnn-c  L2AP for Tanimoto using only Cosine pruning on same threshold",
"  l2ap-m, tapnn-m  L2AP for Tanimoto using MMJoin Tanimoto threshold",
" ",
" utility modes:",
"  info    Get information about the sparse matrix in input-file (output-file ignored).",
"  testeq  Test whether matrix in input-file is the same as that in output-file.",
"          Differences will be printed out to stdout.",
"  io      Transform sparse matrix in input file and write to output-file in",
"          specified format. Scale and Norm parameters can also be invoked.",
" ",
" <input-file> should be in CSR, CLUTO, IJV, AllPairs binary, or binary CSR format.",
" If no <output-file> is given, output will be printed to stdout. In this case, -fmtWrite ",
" must be either IJV (default) or non-binary if -nim is invoked.",
" ",
" Input is assumed to have unit-length rows for Cosine APSS. Otherwise, use the -norm and",
" optionally the -scale parameters to pre-process your input before similarity search.",
" ",
" Options",
"  -t=float",
"     Specifies the similarity threshold used for the search. Should be in (0,1].",
"     Default value is 0.5.",
" ",
"  -sim=string",
"     Which similarity function to use (cos or tan).",
"     Default value is cos.",
" ",
"  -nim",
"     Store neighbors in memory. For some matrices this may produce faster results,",
"     but may require memory many times the input size. See README file for details.",
" ",
"  -scale",
"     Scale the input data by IDF.",
" ",
"  -norm=int",
"     Normalize the matrix rows using the l1 (norm=1) or l2 (norm=2) norm.",
" ",
"  -pr,-pc",
"     Prune rows/cols from the input matrix that are too short/long.",
" ",
"  -prmin=int,-prmax=int,-pcmin=int,-pcmax=int",
"     Minimum/maximum row/column length (nnzs) when pruning (only used with -pr/-pc).",
" ",
"  -fmtRead=string",
"     What format is the dataset stored in: clu, csr, met, ijv, binr, binc, bijv, sbin.",
"     See README for format definitions.",
"     Default value is 0 (detect from extension).",
" ",
"  -readZidx",
"     Column ids start with 0 instead of 1. Pertains to clu, csr, met, and ijv formats only.",
" ",
"  -readVals=int",
"     Read values from file. Pertains to io mode and clu, csr, met, and ijv formats only.",
"     Default value is 1.",
" ",
"  -fmtWrite=string",
"     What format should the output file be written in. See -fmtRead for values.",
"     Default value is ijv.",
" ",
"  -writeZidx",
"     Column ids start with 0 instead of 1. Pertains to clu, csr, met, and ijv formats only.",
" ",
"  -writeVals=int",
"     Write values to file. Pertains to io mode and clu, csr, met, and ijv formats only.",
"     Default value is 1.",
" ",
"  -compactCols, -compactRows",
"     Remove empty cols/rows from the matrix.",
" ",
"  -seed=int",
"     Seed for the random number generator.",
"     Default value is time(NULL).",
" ",
"  -fldelta=int",
"     Float delta used when testing equality of real numbers. (testeq mode only)",
"     Default value is 1e-4.",
" ",
"  -verb=int",
"     Specifies the level of debugging information to be displayed:",
"         0 = NONE, 1 = INFO",
"     Default value is 0 (NONE).",
" ",
"  -v",
"     Prints version information.",
" ",
"  -help, -h",
"     Prints this message.",
""
};


const gk_StringMap_t mode_options[] = {
  {"ij",                MODE_IDXJOIN},
  {"ap",                MODE_AP},
  {"ap2",               MODE_AP2},
  {"l2ap",              MODE_L2AP},
  {"tapnn",             MODE_L2AP},
  {"l2ap-t",            MODE_L2AP_T2},
  {"l2ap-c",            MODE_L2AP_CT},
  {"tapnn-c",           MODE_L2AP_CT},
  {"l2ap-m",            MODE_L2AP_MT},
  {"tapnn-m",           MODE_L2AP_MT},
  {"mmj",               MODE_MMJOIN},
  {"kj",                MODE_MKJOIN},
  {"mkj",               MODE_MKJOIN},
  {"mkj2",              MODE_MKJOIN2},
  {"testeq",            MODE_TESTEQUAL},
  {"io",                MODE_IO},
  {"info",              MODE_INFO},
  {NULL,                0}
};


const gk_StringMap_t sim_options[] = {
  {"cos",               DA_SIM_COS},
  {"jac",               DA_SIM_JAC},
  {"min",               DA_SIM_MIN},
  {"amin",              DA_SIM_AMIN},
  {"euc",               DA_SIM_EUC},
  {"dice",              DA_SIM_DICE},
  {"over",              DA_SIM_OVER},
  {"man",               DA_SIM_MAN},
  {"tan",               DA_SIM_TAN},
  {NULL,                 0}
};

const gk_StringMap_t fmt_options[] = {
  {"clu",               DA_FMT_CLUTO},
  {"csr",               DA_FMT_CSR},
  {"met",               DA_FMT_METIS},
  {"ijv",               DA_FMT_IJV},
  {"binr",              DA_FMT_BINROW},
  {"binc",              DA_FMT_BINCOL},
  {"bijv",              DA_FMT_BIJV},
  {"smat",              DA_FMT_SMAT},
  {"sbin",              DA_FMT_BINAP},
  {"satubin",           DA_FMT_BINAP},
  {"sbinb",             DA_FMT_BINAPB},
  {"satubinb",          DA_FMT_BINAPB},
  {NULL,                 0}
};


const gk_StringMap_t scale_options[] = {
  {"maxtf",            DA_SCALE_MAXTF},
  {"maxtf2",           DA_SCALE_MAXTF2},
  {"sqrt",             DA_SCALE_SQRT},
  {"pow25",            DA_SCALE_POW25},
  {"pow65",            DA_SCALE_POW65},
  {"pow75",            DA_SCALE_POW75},
  {"pow85",            DA_SCALE_POW85},
  {"log",              DA_SCALE_LOG},
  {"idf",              DA_SCALE_IDF},
  {"idf2",             DA_SCALE_IDF2},
  {NULL,               0}
};


/*************************************************************************/
/*! This is the entry point of the command-line argument parser          */
/*************************************************************************/
void cmdline_parse(params_t *params, int argc, char *argv[])
{
	gk_idx_t i;
	int32_t c, option_index;


	/* initialize the params data structure */
	memset(params, 0, sizeof(params_t));

	params->verbosity    = 1;
	params->seed         = -1;
	params->scale        = 0;
	params->norm         = 0;
	params->prunerows    = 0;
	params->prminlen     = -1;
	params->prmaxlen     = -1;
	params->prunecols    = 0;
	params->pcminlen     = -1;
	params->pcmaxlen     = -1;
	params->compactcols  = 0;
	params->compactrows  = 0;
	params->nim          = 0;
	params->sim          = DA_SIM_COS;
	params->nqrows       = 25000;
	params->ndrows       = 100000;
	params->ninnzs       = 1e+6;

	params->simT         = 0.5;
	params->fldelta      = 1e-4;

	params->iFile        = NULL;
	params->fmtRead      = 0;
	params->readVals     = 1;
	params->readNum      = 1;
	params->oFile        = NULL;
	params->fmtWrite     = 0;
	params->writeVals    = 1;
	params->writeNum     = 1;

	params->dataset      = NULL;
	params->filename     = da_cmalloc(1024, "cmdline_parse: filename");
	params->neighcache   = NULL;
	params->neighbors    = NULL;
	params->docs         = NULL;
	params->rperm        = NULL;
	params->cperm        = NULL;
	params->fpout        = NULL;

	/* some defaults */
	params->nSimPairs       = 0;
	params->nCandidates     = 0;
	params->nDotProducts    = 0;
	params->nPruneMinsize   = 0;
	params->nPruneDotP      = 0;
	params->nPruneDotP2     = 0;
	params->nPrunePscore    = 0;
	params->nPruneLength    = 0;
	params->nPruneLength2   = 0;
	params->indexSize       = 0;
	params->nghnnz          = 0;
	params->nghinccnt       = 0;


	/* timers */
	params->timer_1      = 0.0;
	params->timer_2      = 0.0;
	params->timer_3      = 0.0;
	params->timer_4      = 0.0;
	params->timer_5      = 0.0;
	params->timer_6      = 0.0;
	params->timer_7      = 0.0;
	params->timer_8      = 0.0;
	params->timer_9      = 0.0;
	params->timer_global = 0.0;

	/* Parse the command line arguments  */
	while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
		switch (c) {

		case CMD_SIMTHRESHOLD:
			if (gk_optarg) {
				if ((params->simT = atof(gk_optarg)) < 0 || params->simT > 1)
					gk_errexit(SIGERR, "The -simT must be in (0,1].\n");
			}
			break;

        case CMD_SIM:
            if (gk_optarg) {
                if ((params->sim = da_getStringID(sim_options, gk_optarg)) == -1)
                    gk_errexit(SIGERR, "Invalid -sim. Options are: cos or tan.\n");
            }
            break;

		case CMD_FLDELTA:
			if (gk_optarg) {
				if ((params->fldelta = atof(gk_optarg)) <= 0)
					gk_errexit(SIGERR, "The -fldelta must be greater than 0.\n");
			}
			break;


		case CMD_VERBOSITY:
			if (gk_optarg) {
				if ((params->verbosity = atoi(gk_optarg)) < 0)
					gk_errexit(SIGERR, "The -verbosity must be non-negative.\n");
			}
			break;

		case CMD_FMT_READ:
			if (gk_optarg) {
				if ((params->fmtRead = da_getStringID(fmt_options, gk_optarg)) == -1)
					gk_errexit(SIGERR, "Invalid -fmtRead. Options are: clu, csr, met, binr, and binc.\n");
			}
			break;

		case CMD_READ_VALS:
			if (gk_optarg) {
				if ((params->readVals = atoi(gk_optarg)) < 0 || params->readVals > 1)
					gk_errexit(SIGERR, "Invalid -readVals. Must be 0 or 1.\n");
			}
			break;

		case CMD_FMT_READ_NUM:
			params->readNum = 0;
			break;


		case CMD_FMT_WRITE:
			if (gk_optarg) {
				if ((params->fmtWrite = da_getStringID(fmt_options, gk_optarg)) == -1)
					gk_errexit(SIGERR, "Invalid -fmtWrite. Options are: clu, csr, met, binr, and binc.\n");
			}
			break;

		case CMD_WRITE_VALS:
			if (gk_optarg) {
				if ((params->writeVals = atoi(gk_optarg)) < 0 || params->writeVals > 1)
					gk_errexit(SIGERR, "Invalid -writeVals. Must be 0 or 1.\n");
			}
			break;

		case CMD_FMT_WRITE_NUM:
			params->writeNum = 0;
			break;

		case CMD_COMPACT_COLS:
			params->compactcols = 1;
			break;

		case CMD_COMPACT_ROWS:
			params->compactrows = 1;
			break;

		case CMD_SEED:
			if (gk_optarg) {
				if ((params->seed = atoi(gk_optarg)) < 0)
					gk_errexit(SIGERR, "The -seed must be positive.\n");
			}
			break;

		case CMD_SCALE:
			params->scale = 1;
			break;

		case CMD_NIM:
			params->nim = 1;
			break;

		case CMD_PR:
			params->prunerows = 1;
			break;

		case CMD_PRMINLEN:
			if (gk_optarg) {
				if ((params->prminlen = atoi(gk_optarg)) < 1)
					gk_errexit(SIGERR, "The -prmin must be non-negative.\n");
			}
			break;

		case CMD_PRMAXLEN:
			if (gk_optarg) {
				if ((params->prmaxlen = atoi(gk_optarg)) < 1)
					gk_errexit(SIGERR, "The -prmax must be non-negative.\n");
			}
			break;

		case CMD_PC:
			params->prunecols = 1;
			break;

		case CMD_PCMINLEN:
			if (gk_optarg) {
				if ((params->pcminlen = atoi(gk_optarg)) < 1)
					gk_errexit(SIGERR, "The -pcmin must be non-negative.\n");
			}
			break;

		case CMD_PCMAXLEN:
			if (gk_optarg) {
				if ((params->pcmaxlen = atoi(gk_optarg)) < 1)
					gk_errexit(SIGERR, "The -pcmax must be non-negative.\n");
			}
			break;

		case CMD_NORM:
			if (gk_optarg) {
				if ((params->norm = atoi(gk_optarg)) != 1 && params->norm != 2)
					gk_errexit(SIGERR, "The -norm must be 1 or 2.\n");
			}
			break;

		case CMD_VERSION:
			printf("%s (%d.%d.%d), vInfo: [%s]\n", argv[0], VER_MAJOR, VER_MINOR,
						VER_SUBMINOR, VER_COMMENT);
			exit(EXIT_SUCCESS);
			break;

		case CMD_HELP:
			for (i=0; strlen(helpstr[i]) > 0; i++)
				printf("%s\n", helpstr[i]);
			exit(EXIT_SUCCESS);
			break;

		default:
			printf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
			exit(EXIT_FAILURE);
		}
	}

	/* Get the operation to be performed */
	if(argc > gk_optind)
		params->mode      = da_getStringID(mode_options, argv[gk_optind]);
	if (params->mode == -1)
		gk_errexit(SIGERR, "Invalid mode %s.\n", argv[gk_optind]);
	gk_optind++;

	if(argc > gk_optind)
		params->iFile     = gk_strdup(argv[gk_optind++]);
	else
		gk_errexit(SIGERR, "Input file missing.");

	if(!gk_fexists(params->iFile))
		gk_errexit(SIGERR, "Invalid input file %s!\n", params->iFile);

	params->fmtRead = da_getFileFormat(params->iFile, params->fmtRead);

	if(argc > gk_optind){
		/* output file passed in */
		if(params->mode != MODE_IO && params->mode != MODE_TESTEQUAL){
			if(params->nim || da_getFileFormat(argv[gk_optind], -1) == DA_FMT_IJV
			        || gk_strcasecmp(argv[gk_optind], "none") == 0
                    || gk_strcasecmp(argv[gk_optind], "null") == 0){
				params->oFile     = gk_strdup(argv[gk_optind++]);
			} else {
				/* file neighbor output and specified format is not ijv */
				params->oFile     = da_cmalloc(strlen(argv[gk_optind])+9, "cmdline: params->oFile");
				sprintf(params->oFile, "%s.ijv", argv[gk_optind++]);
				if(params->fmtWrite > 0 && params->fmtWrite != DA_FMT_IJV)
					fprintf(stderr, "Warning: Output format reset to IJV. You must invoke -nim to use other output formats.\n");
			}
			if(strcmp(params->iFile, params->oFile) == 0)
				gk_errexit(SIGERR, "The input file and the output file cannot be the same.");
		} else
			params->oFile     = gk_strdup(argv[gk_optind++]);

		params->fmtWrite = da_getFileFormat(params->oFile, params->fmtWrite);

	}

	if(!params->oFile){
		/* no output file name */
		if(params->mode == MODE_TESTEQUAL)
			gk_errexit(SIGERR, "Output file required for mode %s!\n", da_getStringKey(mode_options, params->mode));

		// turn off visibility as output is to stdout
		params->verbosity = -1;

		if(params->mode != MODE_INFO && params->mode != MODE_TESTEQUAL) {
			if(da_isFmtBinary(params->fmtWrite))
				gk_errexit(SIGERR, "Output file must be provided for binary output.");

			if(params->fmtWrite == 0)
				params->fmtWrite = DA_FMT_IJV;

			if( params->mode != MODE_IO && !params->nim && params->fmtWrite != DA_FMT_IJV)
				gk_errexit(SIGERR, "Output format must be IJV unless invoking -nim.");

			params->fpout = stdout;
		}

	}

	if(params->seed < 0)
		params->seed = time(NULL);
	gk_randinit(params->seed);

	if(params->prunerows && params->prminlen == -1 && params->prmaxlen == -1)
		gk_errexit(SIGERR, "You must specify -prmin and/or -prmax when invoking -pr!\n");

	if(params->prunecols && params->pcminlen == -1 && params->pcmaxlen == -1)
		gk_errexit(SIGERR, "You must specify -pcmin and/or -pcmax when invoking -pc!\n");


	/* print the command line */
	if(params->verbosity > 0){
		for (i=0; i<argc; i++)
			printf("%s ", argv[i]);
		printf("\n");
	}
}


