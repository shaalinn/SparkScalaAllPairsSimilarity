/*!
 \file  cmdline.c
 \brief This file contains functions for parsing command-line arguments

 \author David C. Anastasiu
 */
#include "includes.h"


/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct da_option long_options[] = {
    {"k",                 1,      0,      CMD_K},
    {"alpha",             1,      0,      CMD_ALPHA},
    {"alpha",             1,      0,      CMD_ALPHA},
    {"mu",                1,      0,      CMD_MU},
    {"enh",               1,      0,      CMD_NENHANCE},
    {"step",              1,      0,      CMD_STEP},
    {"nbl",               1,      0,      CMD_NBL},
	{"nim",               0,      0,      CMD_NIM},
    {"permcols",          1,      0,      CMD_PERM_COLS},
    {"fldelta",           1,      0,      CMD_FLDELTA},
    {"verb",              1,      0,      CMD_VERBOSITY},
    {"version",           0,      0,      CMD_VERSION},
	{"v",                 1,      0,      CMD_VERIFY},
    {"pr",                0,      0,      CMD_PR},
	{"prmin",             1,      0,      CMD_PRMINLEN},
	{"prmax",             1,      0,      CMD_PRMAXLEN},
    {"pc",                0,      0,      CMD_PC},
    {"cr",                0,      0,      CMD_CR},
	{"pcmin",             1,      0,      CMD_PCMINLEN},
	{"pcmax",             1,      0,      CMD_PCMAXLEN},
	{"fmtRead",           1,      0,      CMD_FMT_READ},
	{"readZidx",          0,      0,      CMD_FMT_READ_NUM},
	{"readVals",          1,      0,      CMD_READ_VALS},
	{"fmtWrite",          1,      0,      CMD_FMT_WRITE},
	{"writeZidx",         0,      0,      CMD_FMT_WRITE_NUM},
	{"writeVals",         1,      0,      CMD_WRITE_VALS},
    {"compactCols",       0,      0,      CMD_COMPACT_COLS},
    {"sortCols",          0,      0,      CMD_SORT_COLS},
    {"compactRows",       0,      0,      CMD_COMPACT_ROWS},
    {"cmpi",              0,      0,      CMD_CMPIND},
    {"cmpv",              0,      0,      CMD_CMPVAL},
    {"stats",             0,      0,      CMD_STATS},
	{"seed",              1,      0,      CMD_SEED},
	{"scale",             0,      0,      CMD_SCALE},
	{"norm",              1,      0,      CMD_NORM},
    {"nqrows",            1,      0,      CMD_NQROWS},
    {"ndrows",            1,      0,      CMD_NDROWS},
    {"nthreads",          1,      0,      CMD_NTHREADS},
	{"help",              0,      0,      CMD_HELP},
	{"h",                 0,      0,      CMD_HELP},
	{0,                   0,      0,      0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
PROGRAM_NAME " - Compute the exact or approximate Cosine K-Nearest Neighbor graph for a set of sparse vectors.",
" ",
"Usage: " PROGRAM_NAME " [options] mode input-file [output-file]",
" ",
" mode:",
"  kij     Build graph using IdxJoin (full sparse dot-product with lesser id docs)",
"  pkij    A multi-threaded version of kij.",
"  kl2ap   Build graph using L2AP (AllPairs Similarity Search with L2-Norm based pruning)",
"  msc     Maxscore Information Retrieval solution",
"  bmm     Block-Max Maxscore with variable block size",
"  bmmc    Block-Max Maxscore with variable block size and compression",
"  gf      Greedy Filtering approximate solution",
"  l2knn   Efficient K-Nearest Neighbor Graph using L2-Norm pruning bounds",
"  l2knn-a Approximate l2knn solution",
" ",
" utility modes:",
"  info    Get information about the sparse matrix in input-file (output-file ignored).",
"  testeq  Test whether matrix in input-file is the same as that in output-file.",
"          Differences will be printed out to stdout.",
"  io      Transform sparse matrix in input file and write to output-file in",
"          specified format. Scale and Norm parameters can also be invoked.",
"  recall  Compute recall of a knng solution given true values.",
" ",
" <input-file> should be in CSR, CLUTO, IJV, AllPairs binary, or binary CSR format.",
" Specify stdout for the <output-file> to print results to stdout. In this case, -fmtWrite ",
" must be either IJV (default) or non-binary if -nim is invoked. If no <output-file> is",
" specified, the output will not be saved. K-NNG output will always be sparse vectors, ",
" sorted in decreasing similarity order.",
" ",
" Input is assumed to have unit-length rows when computing cosine similarity. Otherwise, use",
" the -norm and optionally the -scale parameters to pre-process input before similarity search.",
" ",
" Options",
"  -k=int",
"     Number of neighbors to return for each row in the K-Nearest Neighbor Graph.",
"     Default value is 10.",
" ",
"  -alpha=float",
"     Number of neighbors to consider initially, as a multiple of k.",
"     Default value is 2, i.e. will consider 2*k. Must be non-negative.",
" ",
"  -mu=float",
"     Number of neighbors to consider initially. If set, alpha = mu/k.",
"     Default value is alpha*k. Must be >= k.",
" ",
"  -enh=int",
"     Number of iterations for the neighborhood enhancing algo.",
"     Default value is 10.",
" ",
"  -step=float",
"     Step size to decrease simT by for kl2ap mode.",
"     Default value is 0.1. Must be in (0,1].",
" ",
"  -nbl=int",
"     Number of blocks to split the K-Nearest Neighbor Graph processing in mode l2knn.",
"     Default value is 1 (no splitting).",
" ",
"  -v=string",
"     Verification file containing a true K-Nearest Neighbor Graph. Must be in CSR format.",
"     Default value is NULL (no verification).",
" ",
"  -permcols=string",
"     How to permute columns before computing knng (for kij and pkij modes): none, l2ap, or l2knn.",
"     Default: none.",
" ",
"  -cr",
"     Use checkpoint-restart for methods that have the technology enabled.",
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
"  -cmpi,-cmpv",
"     Do not compare indices/values when comparing matrices in mode testeq. ",
"     Defaults to compare both indices and values.",
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
"  -sortCols",
"     Sort column indices in matrix.",
" ",
"  -stats",
"     Display additional statistics for the matrix (applies to mode 'info' only).",
" ",
"  -seed=int",
"     Seed for the random number generator.",
"     Default value is time(NULL).",
" ",
"  -fldelta=int",
"     Float delta used when testing equality of real numbers. (testeq mode only)",
"     Default value is 1e-4.",
" ",
"  -nthreads=int",
"     The number of threads per block to be used for computing the neighbors.",
"     Default value is 1.",
" ",
"  -verb=int",
"     Specifies the level of debugging information to be displayed:",
"         0 = NONE, 1 = INFO",
"     Default value is 0 (NONE).",
" ",
"  -version",
"     Prints version information.",
" ",
"  -help, -h",
"     Prints this message.",
""
};


const da_StringMap_t mode_options[] = {
  {"kij",               MODE_KIDXJOIN},
  {"kl2ap",             MODE_KL2AP},
  {"l2knn",             MODE_L2KNN},
  {"l2knn-a",           MODE_L2KNN_A},
  {"approx",            MODE_L2KNN_A},
  {"pkij",              MODE_PKIDXJOIN},
  {"msc",               MODE_MSC},
  {"bmm",               MODE_BMM},
  {"bmmc",              MODE_BMMC},
  {"gf",                MODE_GF},
  {"recall",            MODE_RECALL},
  {"testeq",            MODE_TESTEQUAL},
  {"io",                MODE_IO},
  {"info",              MODE_INFO},
  {NULL,                 0}
};


const da_StringMap_t sim_options[] = {
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


const da_StringMap_t perm_options[] = {
  {"none",               PERM_NONE},
  {"l2ap",               PERM_L2AP},
  {"l2knn",               PERM_L2KNN},
  {NULL,                 0}
};

const da_StringMap_t fmt_options[] = {
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


const da_StringMap_t scale_options[] = {
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
	idx_t i;
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
	params->sortcols     = 0;
    params->cmpind       = 1;
    params->cmpval       = 1;
	params->stats        = 0;
	params->nim          = 0;
    params->nqrows       = 5000;
    params->ndrows       = 1000;
    params->permcol      = PERM_NONE;
    params->nthreads     = 1;
    params->nbl          = 1;
    params->step         = 0.1;
    params->cr           = 0;

	params->sim          = DA_SIM_COS;
	params->simT         = 0.5;
	params->fldelta      = 1e-4;
	params->k            = 10;
    params->alpha        = 2;
    params->mu           = 0;
	params->pmerge       = PMERGE_TWOS;
	params->nenhance     = 2;

	params->iFile        = NULL;
	params->fmtRead      = 0;
	params->readVals     = 1;
	params->readNum      = 1;
	params->oFile        = NULL;
	params->fmtWrite     = 0;
	params->writeVals    = 1;
	params->writeNum     = 1;
    params->vFile        = NULL;

	params->dataset      = NULL;
	params->filename     = da_cmalloc(1024, "cmdline_parse: filename");
    params->crfname      = da_cmalloc(1024, "cmdline_parse: crfname");
    params->crfname2     = da_cmalloc(1024, "cmdline_parse: crfname");
	params->neighcache   = NULL;
	params->neighbors    = NULL;
	params->querydocs    = NULL;
	params->knng         = NULL;
	params->docs         = NULL;
	params->rperm        = NULL;
	params->cperm        = NULL;
    params->rsizes       = NULL;
    params->csizes       = NULL;
    params->rwgts        = NULL;
    params->cwgts        = NULL;
	params->fpout        = NULL;

	/* some defaults */
	params->nSimPairs       = 0;
	params->nCandidates     = 0;
	params->nDotProducts    = 0;
    params->nDotProducts1   = 0;
    params->nDotProducts2   = 0;
	params->nPruneMinsize   = 0;
	params->nPruneDotP      = 0;
	params->nPruneDotP2     = 0;
	params->nPrunePscore    = 0;
	params->nPruneLength    = 0;
	params->nPruneLength2   = 0;
	params->accumPriorPrune = 0;
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
	while ((c = da_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
		switch (c) {

        case CMD_PERM_COLS:
            if (da_optarg) {
                if ((params->permcol = da_getStringID(perm_options, da_optarg)) == -1)
                    da_errexit("Invalid -permcols. Options are: none, l2ap, l2knn.\n");
            }
            break;

        case CMD_K:
            if (da_optarg) {
                if ((params->k = atoi(da_optarg)) < 1)
                    da_errexit("Invalid -k. Must be greater than 0.\n");
            }
            break;

        case CMD_NBL:
            if (da_optarg) {
                if ((params->nbl = atoi(da_optarg)) < 1)
                    da_errexit("Invalid -nbl. Must be greater than 0.\n");
            }
            break;

        case CMD_NQROWS:
          if (da_optarg) {
            if ((params->nqrows = atoi(da_optarg)) < 0)
              da_errexit("The -nqrows value must be greater than 0.\n");
          }
          break;

        case CMD_NDROWS:
          if (da_optarg) {
            if ((params->ndrows = atoi(da_optarg)) < 0)
              da_errexit("The -ndrows value must be greater than 0.\n");
          }
          break;

        case CMD_ALPHA:
            if (da_optarg) {
                if ((params->alpha = atof(da_optarg)) < 0)
                    da_errexit("The -alpha value must be non-negative.\n");
            }
            break;

        case CMD_MU:
            if (da_optarg) {
                if ((params->mu = atof(da_optarg)) < 0)
                    da_errexit("The -mu value must be non-negative.\n");
            }
            break;

        case CMD_NENHANCE:
            if (da_optarg) {
                if ((params->nenhance = atoi(da_optarg)) < 0)
                    da_errexit("Invalid -enh. Must be non-negative.\n");
            }
            break;

        case CMD_FLDELTA:
            if (da_optarg) {
                if ((params->fldelta = atof(da_optarg)) <= 0)
                    da_errexit("The -fldelta value must be greater than 0.\n");
            }
            break;

        case CMD_STEP:
            if (da_optarg) {
                if ((params->step = atof(da_optarg)) <= 0 || params->step > 1)
                    da_errexit("The -step value must be in (0,1].\n");
            }
            break;


		case CMD_VERBOSITY:
			if (da_optarg) {
				if ((params->verbosity = atoi(da_optarg)) < 0)
					da_errexit("The -verbosity value must be non-negative.\n");
			}
			break;

		case CMD_FMT_READ:
			if (da_optarg) {
				if ((params->fmtRead = da_getStringID(fmt_options, da_optarg)) == -1)
					da_errexit("Invalid -fmtRead. Options are: clu, csr, met, binr, and binc.\n");
			}
			break;

		case CMD_READ_VALS:
			if (da_optarg) {
				if ((params->readVals = atoi(da_optarg)) < 0 || params->readVals > 1)
					da_errexit("Invalid -readVals. Must be 0 or 1.\n");
			}
			break;

		case CMD_FMT_READ_NUM:
			params->readNum = 0;
			break;


		case CMD_FMT_WRITE:
			if (da_optarg) {
				if ((params->fmtWrite = da_getStringID(fmt_options, da_optarg)) == -1)
					da_errexit("Invalid -fmtWrite. Options are: clu, csr, met, binr, and binc.\n");
			}
			break;

		case CMD_WRITE_VALS:
			if (da_optarg) {
				if ((params->writeVals = atoi(da_optarg)) < 0 || params->writeVals > 1)
					da_errexit("Invalid -writeVals. Must be 0 or 1.\n");
			}
			break;

		case CMD_FMT_WRITE_NUM:
			params->writeNum = 0;
			break;

        case CMD_COMPACT_COLS:
            params->compactcols = 1;
            break;

        case CMD_SORT_COLS:
            params->sortcols = 1;
            break;

        case CMD_COMPACT_ROWS:
            params->compactrows = 1;
            break;

        case CMD_STATS:
            params->stats = 1;
            break;

		case CMD_SEED:
			if (da_optarg) {
				if ((params->seed = atoi(da_optarg)) < 0)
					da_errexit("The -seed value must be positive.\n");
			}
			break;

		case CMD_SCALE:
			params->scale = 1;
			break;

        case CMD_CR:
            params->cr = 1;
            break;

		case CMD_NIM:
			params->nim = 1;
			break;

		case CMD_PR:
			params->prunerows = 1;
			break;

		case CMD_PRMINLEN:
			if (da_optarg) {
				if ((params->prminlen = atoi(da_optarg)) < 1)
					da_errexit("The -prmin value must be non-negative.\n");
			}
			break;

		case CMD_PRMAXLEN:
			if (da_optarg) {
				if ((params->prmaxlen = atoi(da_optarg)) < 1)
					da_errexit("The -prmax value must be non-negative.\n");
			}
			break;

        case CMD_PC:
            params->prunecols = 1;
            break;

        case CMD_CMPIND:
            params->cmpind = 0;
            break;

        case CMD_CMPVAL:
            params->cmpval = 0;
            break;

		case CMD_PCMINLEN:
			if (da_optarg) {
				if ((params->pcminlen = atoi(da_optarg)) < 1)
					da_errexit("The -pcmin value must be non-negative.\n");
			}
			break;

		case CMD_PCMAXLEN:
			if (da_optarg) {
				if ((params->pcmaxlen = atoi(da_optarg)) < 1)
					da_errexit("The -pcmax value must be non-negative.\n");
			}
			break;

		case CMD_NORM:
			if (da_optarg) {
				if ((params->norm = atoi(da_optarg)) != 1 && params->norm != 2)
					da_errexit("The -norm value must be 1 or 2.\n");
			}
			break;

	      case CMD_NTHREADS:
	        if (da_optarg) {
	          if ((params->nthreads = atoi(da_optarg)) < 1)
	            da_errexit("The -nthreads must be greater than 1.\n");
	        }
	        break;

		case CMD_VERSION:
			printf("%s (%d.%d.%d), vInfo: [%s]\n", argv[0], VER_MAJOR, VER_MINOR,
						VER_SUBMINOR, VER_COMMENT);
			exit(EXIT_SUCCESS);
			break;

        case CMD_VERIFY:
            params->vFile = da_strdup(da_optarg);
            if(!da_fexists(params->vFile))
                da_errexit("The -v parameter requires a valid verification file. %s is not a file.\n", params->vFile);
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
	if(argc > da_optind)
		params->mode      = da_getStringID(mode_options, argv[da_optind]);
	if (params->mode == -1)
		da_errexit("Invalid mode %s.\n", argv[da_optind]);
	da_optind++;

	if(argc > da_optind)
		params->iFile     = da_strdup(argv[da_optind++]);
	else
		da_errexit("Input file missing.");

	if(!da_fexists(params->iFile))
		da_errexit("Invalid input file %s!\n", params->iFile);

	params->fmtRead = da_getFileFormat(params->iFile, params->fmtRead);

	if(argc > da_optind){		/* output file passed in */
        params->oFile     = da_strdup(argv[da_optind++]);
        if(strcmp(params->iFile, params->oFile) == 0)
            da_errexit("The input file and the output file cannot be the same.");

		if(da_strcasecmp(params->oFile, "stdout") || da_strcasecmp(params->oFile, "stderr") ||
            da_strcasecmp(params->oFile, "out") || da_strcasecmp(params->oFile, "err")){

            // turn off visibility as output is to stdout or stderr
            params->verbosity = -1;

            if(params->mode != MODE_INFO && params->mode != MODE_TESTEQUAL && params->mode != MODE_RECALL) {
                if(da_isFmtBinary(params->fmtWrite))
                    da_errexit("Output file must be provided for binary output.");

                if(params->fmtWrite == 0)
                    params->fmtWrite = DA_FMT_IJV;

                if( params->mode != MODE_IO && !params->nim && params->fmtWrite != DA_FMT_IJV)
                    da_errexit("Output format must be IJV unless invoking -nim.");

                params->fpout = da_strcasecmp(params->oFile, "stdout") ||
                        da_strcasecmp(params->oFile, "out") ? stdout : stderr;
            }

        } else if(params->mode != MODE_INFO && params->mode != MODE_TESTEQUAL && params->mode != MODE_RECALL) {
            params->fmtWrite = da_getFileFormat(params->oFile, params->fmtWrite);
            params->fpout = da_fopen(params->oFile, "w", "cmdline: oFile fpout");
        }
	}


	if(!params->oFile && params->mode == MODE_TESTEQUAL)
        da_errexit("Output file required for mode %s!\n", da_getStringKey(mode_options, params->mode));

	if(params->mu > 0){
	    if(params->mu < params->k)
	        da_errexit("The -mu value must be >= -k.");
	    params->alpha = params->mu / params->k;
	} else
	    params->mu = params->alpha * params->k;

	if(params->seed < 0)
		params->seed = time(NULL);
	da_randinit(params->seed);

	if(params->prunerows && params->prminlen == -1 && params->prmaxlen == -1)
		da_errexit("You must specify -prmin and/or -prmax when invoking -pr!\n");

	if(params->prunecols && params->pcminlen == -1 && params->pcmaxlen == -1)
		da_errexit("You must specify -pcmin and/or -pcmax when invoking -pc!\n");

	/* print the command line */
	if(params->verbosity > 0){
		for (i=0; i<argc; i++)
			printf("%s ", argv[i]);
		printf("\n");
	}
}


