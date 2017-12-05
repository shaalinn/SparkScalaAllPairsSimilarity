/*!
 \file  main.c
 \brief This file is the entry point for apss's various components

 \author David C. Anastasiu
 */

#include "includes.h"
#include <string.h>

/** forward defs */
void printDefinedOptions (const params_t * const params);

/***********************************************/
/*! This is the entry point for the program    */
/***********************************************/


long lineCount(char *file){
	FILE *fp;
	fp=fopen(file,"r");
	long count=0;
	char c;
	if(fp==NULL){
		printf("\nCan not read file %s",file);
		return -1;
	}
	for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n') // Increment count if this character is newline
            count = count + 1;
    fclose(fp);
    printf("\ncount=%ld",count);
    return count;
}




int main (int argc, char *argv[])
{
    double timer;
    params_t *params;
    time_t seconds;
	seconds=time(NULL);
	struct timeval now;
	gettimeofday(&now, NULL);

	srand(time(NULL));   // should only be called once
	int r = rand();



    FILE *input,*part;
	char file[512],fname[2][512];
	char inputname[512],opfile[512];
	char ch;
	int file_count=0;
	sprintf(inputname,"/home/hadoop/SparkScalaAllPairsSimilarity/pl2ap/generated/input%ld-%d",now.tv_usec,r);
	sprintf(opfile,"/home/hadoop/SparkScalaAllPairsSimilarity/pl2ap/generated/output%ld-%d.ijv",now.tv_usec,r);
	input=fopen(inputname,"w+");
	while(1){
		char ch=fgetc(stdin);
		
		if(ch==EOF){
			break;		
		}
		//printf("%c",ch);	
		fputc(ch,input);
	}
	fclose(input);
	input=fopen(inputname,"r");


		
	//printf("%ld %ld",now.tv_usec,r);
	sprintf(fname[0],"/home/hadoop/SparkScalaAllPairsSimilarity/pl2ap/generated/part-test%ld-%d.csr",now.tv_usec,r);
	sprintf(fname[1],"/home/hadoop/SparkScalaAllPairsSimilarity/pl2ap/generated/part-test%ld-%d(1).csr",now.tv_usec,r);
	if(input == NULL)
    {
        printf("File Not Found\n");
        exit(0);
    }
    else
    {
        while(1)
        {
            ch = fgetc(input);
            if(ch == EOF)   
            {
                break;
            }
            if(ch == ' ')
                file_count++;
            else if(ch == '\t')
                file_count++;
            else if(ch == '\n')
                file_count++;
        }
    }
    fclose(input);
    int i,j;
    if(file_count<1){
    	printf("Invalid input file");
    	exit(0);
    }
    printf("file count = %d\n",file_count);
    char **filename=(char**)malloc(sizeof(char*)*file_count);
    long **index = (long**)malloc(sizeof(long*)*file_count);
    
    input=fopen(inputname,"r");
    
    for(i=0;i<file_count;i++){
    	int j=0;
    	filename[i]=(char*)malloc(sizeof(char)*512);
    	while(1){
    		ch=fgetc(input);
    		if(ch == ' ' || ch=='\n' || ch == EOF){
    			filename[i][j++]='\0';
    			break;
    		}
    		filename[i][j++]=ch;
		
    	}
    }
    fclose(input);
    
//    FILE *temp;
//	temp=fopen("/home/shalin/Downloads/pl2ap/generated/input","w+");

    for(i=0;i<file_count;i++){
    	puts(filename[i]);
    	
    }
    for(i=0;i<file_count;i++)
    {
    	long count=lineCount(filename[i]);
    	if(count!=-1){
    		index[i]=(long*)malloc(sizeof(long)*(count+1));
    	}else{
    		index[i]=NULL;
    	}
    }
    FILE *fp, *output;
    char * line = NULL;
    char *ptr;
    char threshold[20];
    size_t len = 0;
    ssize_t read;
    int linenum;
    long k;
    for(i=0;i<file_count;i++)
    {
    	fp = fopen(filename[i], "r");
    	output = fopen(fname[i],"w+");
	
    	if (fp == NULL){
        	printf("\nCan not read file %s",file);
        	continue;
    	}
    	//printf("\nArray %d: ",i);
    	j=1;
    	while ((read = getline(&line, &len, fp)) != -1) {
    		linenum=strtol(line, &ptr, 10);
        	index[i][j]=linenum;
		fputs(ptr,output);
        	j++;
    	}

    	fclose(fp);
    	fclose(output);
    }
    fp=fopen("/home/hadoop/efs/storage_temp/cluHeaders.txt","r");
    if (fp == NULL){
        	printf("\nCan not read file %s",file);
        	exit(0);
    	}
    	//printf("\nArray %d: ",i);
    	j=1;
    	while ((read = getline(&line, &len, fp)) != -1) {
    		strcpy(threshold,line);
    	}

    	fclose(fp);
    	


    params = (params_t *) da_malloc (sizeof (params_t), "main: params");
    char *args[]={"./pl2ap","pijnn","-norm","2","-nthreads","1",fname[0],fname[1],opfile,"-eps",threshold};

    printf("before command line");
    cmdline_parse (params, 11, args);

    if (params->verbosity > 0) {
        //printf("********************************************************************************\n");
        //printf ("%s (%d.%d.%d), vInfo: [%s]\n", PROGRAM_NAME, VER_MAJOR,
        //        VER_MINOR, VER_SUBMINOR, VER_COMMENT);
        //printf ("mode: %s, ", da_getStringKey (mode_options, params->mode));
        //printf ("iFile: %s, ", params->iFile);
        //printf ("oFile: %s, ", params->oFile ? params->oFile : "NULL");

        /*switch (params->mode) {

            case MODE_IDXJOIN:
                printf ("eps: %.3f", params->eps);
                break;

            case MODE_PIDXJOIN:
                printf ("eps: %.3f\nnqrows: %d, ndrows: %d, nthreads: %d",
                        params->eps, params->nqrows, params->ndrows,
                        params->nthreads);
                break;

            case MODE_PIJNN:
                printf ("eps: %.3f\nnqrows: %d, ndrows: %d, nthreads: %d",
                        params->eps, params->nqrows, params->ndrows,
                        params->nthreads);
                break;

            case MODE_L2AP:
                printf ("eps: %.3f", params->eps);
                break;

            case MODE_PL2AP:
                printf ("eps: %.3f, nqrows: %d, ndrows: %d, nthreads: %d",
                        params->eps, params->nqrows, params->ndrows,
                        params->nthreads);
                break;

            case MODE_PL2NN:
                printf ("eps: %.3f, nqrows: %d, ndrows: %d, nthreads: %d",
                        params->eps, params->nqrows, params->ndrows,
                        params->nthreads);
                break;

            case MODE_TESTEQUAL:
                printf ("fldelta: %g", params->fldelta);
                break;

            default:
                break;
        }*/
        //printDefinedOptions (params);

        //if (params->run)
            //printf ("\nrun: %s", params->run);
        //printf("\n********************************************************************************\n");
        fflush (stdout);
    }

    da_startwctimer (params->timer_global);

    // read input data
    readInputData (params);

    switch (params->mode) {

        case MODE_IDXJOIN:
            ijFindNeighbors (params);
            break;

        case MODE_PIDXJOIN:
            pijFindNeighbors (params);
            break;

        case MODE_PIJNN:
            pijnnFindNeighbors (params);
            break;

        case MODE_L2AP:
            l2apFindNeighbors (params);
            break;

        case MODE_PL2AP:
            pl2apFindNeighbors (params);
            break;

        case MODE_PL2NN:
            pl2nnFindNeighbors (params);
            break;

        case MODE_TESTEQUAL:
            da_testMatricesEqual (params);
            break;

        case MODE_INFO:
            da_matrixInfo (params);
            break;

        case MODE_IO:
            da_matrixIo (params);
            break;

        case MODE_RECALL:
            da_errexit("Mode recall not yet implemented.");
//            da_testRecall (params);
            break;

        default:
            da_errexit ("Invalid mode.");
            break;
    }


    // similarity search complete.
    if (params->verbosity > 0) {
#ifdef EXTRACOUNTS
        if (params->mode == MODE_L2AP || params->mode == MODE_PL2AP) {
            //printf ("Dynamic index size: %zu\n", params->indexSize);
            //printf ("Nnzs pruned by minsize bound: %zu\n",
            //        params->nPruneMinsize);
            //printf ("Num. candidates pruned by length bound (gen): %zu\n",
            //        params->nPruneLength);
            //printf ("Num. candidates pruned by length bound (ver): %zu\n",
            //        params->nPruneLength2);
            //printf ("Num. candidates pruned by pscore bound: %zu\n",
            //        params->nPrunePscore);
            //printf ("Num. candidates pruned by dotp bound: %zu\n",
            //        params->nPruneDotP);
            //printf ("Num. candidates pruned by positional dotp bound: %zu\n",
            //        params->nPruneDotP2);
        }
#endif

        //printf ("Num. total candidates: %zu\n", params->nCandidates);
        //printf ("Num. dot products: %zu\n", params->nDotProducts);
        //printf ("Num. dot products in cg: %zu\n", params->nDotProducts1);
        //printf ("Num. similar pairs: %zu\n", params->nSimPairs);
        //if (params->nim && params->mode == MODE_L2AP)
            //printf ("Num. times neighbors structure increased in size: "
            //        PRNT_PTRTYPE "\n", params->nghinccnt);
        fflush (stdout);

        da_stopwctimer (params->timer_global);

#ifdef EXTRATIMES
        //printf ("TIMES:\n");
        //da_printTimerLong ("\t Memory alloc time: ",
        //        da_getwctimer (params->timer_5));
        /*if (params->mode != MODE_IDXJOIN && params->mode != MODE_PIDXJOIN) {
            da_printTimerLong ("\t Add to index: ",
                    da_getwctimer (params->timer_4));
            da_printTimerLong ("\t Generate candidates: ",
                    da_getwctimer (params->timer_1));
            da_printTimerLong ("\t Process candidates: ",
                    da_getwctimer (params->timer_2));
            timer = params->timer_3 - params->timer_1 - params->timer_2 -
                    params->timer_4;
            da_printTimerLong ("\t Excess search time: ",
                    da_getwctimer (timer));
        }*/
#endif
        //da_printTimerLong ("\t Similarity search: ",
        //        da_getwctimer (params->timer_3));
        //da_printTimerLong ("\t Total time: ",
        //        da_getwctimer (params->timer_global));

        //printf
        //("********************************************************************************\n");
    }

    freeParams (&params);




	//**************************************   Output processing ***********************************
	fp = fopen(opfile, "r");
    	if (fp == NULL){
        	printf("\nCan not read file %s",file);
 		exit(0);
    	}
    	
    	
    	while ((read = getline(&line, &len, fp)) != -1) {
    		linenum=strtol(line, &ptr, 10);
        	printf("%ld ",index[0][(int)linenum]-1);
		strcpy(line,ptr);
    		linenum=strtol(line, &ptr, 10);
		printf("%ld ",index[1][(int)linenum]-1);
		printf("%s",ptr);
        	
    	}

    	fclose(fp);






    exit (EXIT_SUCCESS);
}


void printDefinedOptions (const params_t * const params)
{

    switch (params->mode)
    {
        case MODE_PL2AP:
#if PL2APRS == PL2APRS_STORE
            printf (", store_rs: true");
#else
            printf (", store_rs: false");
#endif
            printf ("\npruning: min(b1,b3)");
#ifdef PL2AP_MINSZ
#ifdef PL2AP_SUP
            printf (", sz3 (upd)");
#else
            printf (", sz3");
#endif
#endif
            printf (", min(rs1,rs4), l2cg, ps, dp5");
#ifdef PL2AP_DP6
            printf (", dp6");
#endif
            printf (", l2cv");

#if L2APDP == L2APDP_MIX
            printf ("\nhashing: mix");
#elif L2APDP == L2APDP_HASH
            printf ("\nhashing: dense");
#elif L2APDP == L2APDP_MASK
            printf ("\nhashing: mask");
#endif
#if L2APDP == L2APDP_MIX || L2APDP == L2APDP_MASK
            printf (", ht size: %d, ht ub: %d", HTSIZE, HTUBND);
#endif

            break;

        default:
            break;
    }
}


/**
 * Pre-process dataset by pruning, scaling, filtering, shuffling, etc.
 */
void preProcessData (params_t * params)
{
    da_csr_t *docs = params->docs;
    da_csr_t *tmp = NULL;

    /* prune cols with less than pcminlen or more than pcmaxlen values */
    if (params->prunecols) {
        if (params->pcminlen == -1)
            params->pcminlen = 0;
        if (params->pcmaxlen == -1)
            params->pcmaxlen = docs->nrows;
        printf ("pc: %d, pcmin: %d, pcmax: %d\n", params->prunecols,
                params->pcminlen, params->pcmaxlen);
        tmp = da_csr_Prune (docs, DA_COL, params->pcminlen, params->pcmaxlen);
        da_csr_Transfer (tmp, docs);
        da_csr_Free (&tmp);
    }

    /* prune rows with less than prminlen or more than prmaxlen values */
    if (params->prunerows) {
        if (params->prminlen == -1)
            params->prminlen = 0;
        if (params->prmaxlen == -1)
            params->prmaxlen = docs->ncols;
        printf ("pr: %d, prmin: %d, prmax: %d\n", params->prunerows,
                params->prminlen, params->prmaxlen);
        tmp = da_csr_Prune (docs, DA_ROW, params->prminlen, params->prmaxlen);
        da_csr_Transfer (tmp, docs);
        da_csr_Free (&tmp);
    }


    // compact the column space
    if (params->compactcols)
        da_csr_CompactColumns (docs);

    // sort the column space
    if (params->sortcols)
        da_csr_SortIndices (docs, DA_ROW);

    // compact the row space
    if (params->compactrows)
        da_csr_CompactRows (docs);

    /* scale term values */
    if (params->scale) {
        if (params->verbosity > 0)
            printf ("   Scaling input matrix.\n");
        da_csr_Scale (docs, DA_SCALE_IDF);
    }

    /* normalize docs rows */
    if (params->norm > 0)
        da_csr_Normalize (docs, DA_ROWCOL, params->norm);

}

void simSearchSetup (params_t * params)
{
    da_csr_t *docs = params->docs;
    da_csr_t *neighbors = NULL;
    size_t nghnnz;

    // pre-process the data
    preProcessData (params);

    // compact the column space
    da_csr_CompactColumns (docs);

    if (params->verbosity > 0)
        printf ("Docs matrix: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, "
                PRNT_PTRTYPE " nnz\n", docs->nrows, docs->ncols,
                docs->rowptr[docs->nrows]);

#ifndef NO_OUTPUT
    da_startwctimer (params->timer_5);	// memory allocation time
    if (params->nim) {
        // allocate memory for tracking neighbors - O <-- {}
        neighbors = da_csr_Create ();
        neighbors->nrows = neighbors->ncols = docs->nrows;
        nghnnz = params->nghnnz = NINITNEIGHBORS * docs->nrows;
        neighbors->rowptr =
                da_pmalloc (docs->nrows + 1, "apFindNeighbors: neighbors->rowptr");
        neighbors->rowind =
                da_imalloc (nghnnz, "apFindNeighbors: neighbors->rowind");
        neighbors->rowval =
                da_vmalloc (nghnnz, "apFindNeighbors: neighbors->rowval");
        neighbors->rowptr[0] = 0;
        params->neighbors = neighbors;

    } else {
        // set up neighbor output cache
        params->neighcache =
                (da_sims_t *) da_malloc (sizeof (da_sims_t),
                        "l2apFindNeighbors: neighcache");
        params->neighcache->nsims = 0;
        params->neighcache->size = NSIMSBUFFER;
        params->neighcache->sims = da_ssmalloc (NSIMSBUFFER, (da_sim_t) {0, 0, 0.0},
        "l2apFindNeighbors: neighcache->sims");

        /* open output file */
        if (params->oFile && params->mode < MODE_INFO)
            params->fpout = da_fopen (params->oFile, "w", "da_csr_Write: fpout");
    }
    da_stopwctimer (params->timer_5);	// memory allocation time
#endif
}



void simSearchFinalize (params_t * params)
{
#ifndef NO_OUTPUT
    if (params->nim) {
        // inverse permute neighbors
        if (params->rperm) {
            da_inversePermuteSqMatrixLT (params->neighbors, params->rperm,
                    DA_ROW);
        }
        if (params->fpout) {
            // make output symmetric
            if (params->symout)
                da_addHigherNeighbors (params->neighbors);
            if (params->verbosity > 0)
                printf ("Writing neighborhood matrix to %s.\n", params->oFile);
            da_csr_Write (params->neighbors, params->oFile, params->fmtWrite, 1,
                    1);
        }
    } else {
        // write out neighbors buffer
        da_storeSim (params, -1, -1, 0);
    }
#endif
}

/**
 * Read input data
 */
void readInputData (params_t * params)
{
    da_csr_t *docs;
    params->fmtRead = da_getFileFormat (params->iFile, params->fmtRead);
    if (params->fmtRead < 1)
        da_errexit ("Invalid input format.\n");
    docs = da_csr_Read (params->iFile, params->fmtRead, params->readVals,
                    params->readNum);
    assert (docs->rowptr || docs->colptr);
    /* ensure row based structure for sim search */
    if (params->mode < 50 && !docs->rowptr) {
        da_csr_CreateIndex (docs, DA_ROW);
        da_csr_FreeBase (docs, DA_COL);
    }
    params->docs = docs;
    if(params->dFile){
        params->docsdb = da_csr_Read (params->dFile, params->fmtRead, params->readVals,
                params->readNum);
    }
}

/**
 * Test two matrices are equal. Values tested up to params->fldelta precision.
 */
void da_testMatricesEqual (params_t * params)
{
    da_csr_t *docs = NULL, *docs2 = NULL;

    docs = params->docs;

    // test equality of two sparse matrices & print out differences
    docs2 = da_csr_Read (params->oFile,
                    da_getFileFormat (params->oFile, params->fmtWrite),
                    params->readVals, params->readNum);
    printf ("Comparing %s (A[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE
            "]) and " "%s (B[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE
            "]).\n\n", params->iFile, docs->nrows, docs->ncols,
            docs->rowptr[docs->nrows], params->oFile, docs2->nrows,
            docs2->ncols, docs2->rowptr[docs2->nrows]);
    da_csrCompare (docs, docs2, params->fldelta, params->cmpind,
            params->cmpval);
    da_csr_Free (&docs2);
    freeParams (&params);
    exit (EXIT_SUCCESS);
}




/**
 * Display information about a sparse matrix
 */
void da_matrixInfo (params_t * params)
{
    da_csr_t *docs = params->docs;
    size_t i, j, l, x1, x2, min, max;

    // identify format, read matrix, and print out information about it - nrows, ncols, nnz
    printf ("%s: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, " PRNT_PTRTYPE
            " nnzs, %g density, ", params->iFile, docs->nrows, docs->ncols,
            docs->rowptr[docs->nrows],
            docs->rowptr[docs->nrows] / ((double) docs->nrows * docs->ncols));

    da_csr_CompactColumns (docs);
    printf (PRNT_IDXTYPE " non-empty cols.\n", docs->ncols);

    if (params->stats) {
        if (!docs->colptr)
            da_csr_CreateIndex (docs, DA_COL);
        for (x1 = 0, x2 = 0, min = INT_MAX, max = 0, i = 0; i < docs->nrows; ++i) {
            l = docs->rowptr[i + 1] - docs->rowptr[i];
            x1 += l;
            x2 += l * l;
            if (l > max)
                max = l;
            if (l < min)
                min = l;
        }
        printf ("Row nnz stats: min %zu, max %zu mean %.2f, stdev %.2f.\n", min,
                max, (double) x1 / (double) docs->nrows,
                sqrt ((double) x2 * (double) docs->nrows -
                        x1 * x1) / (double) docs->nrows);
        for (x1 = 0, x2 = 0, i = 0; i < docs->ncols; ++i) {
            l = docs->colptr[i + 1] - docs->colptr[i];
            x1 += l;
            x2 += l * l;
            if (l > max)
                max = l;
            if (l < min)
                min = l;
        }
        printf ("Col nnz stats: min %zu, max %zu mean %.2f, stdev %.2f.\n", min,
                max, (double) x1 / (double) docs->ncols,
                sqrt ((double) x2 * (double) docs->ncols -
                        x1 * x1) / (double) docs->ncols);
    }



    printf ("\n");

    freeParams (&params);
    exit (EXIT_SUCCESS);
}


/**
 * Transform input matrix to some other format
 * Data pre-processing invoked before transforming.
 */
void da_matrixIo (params_t * params)
{
    da_csr_t *docs = params->docs;

    if (params->oFile) {
        params->fmtWrite = da_getFileFormat (params->oFile, params->fmtWrite);
        if (params->fmtWrite < 1)
            da_errexit ("Invalid output format.\n");
    }

    if (params->fmtRead == DA_FMT_BINAPB && params->fmtWrite != DA_FMT_BINAPB) {
        // transition from binary to weighted format
        if (!docs->rowval)
            docs->rowval = da_vsmalloc (docs->rowptr[docs->nrows], 1.0, "io: docs->rowptr");
    }

    if (params->verbosity > 0)
        printf ("Transforming %s (A[" PRNT_IDXTYPE "," PRNT_IDXTYPE ","
                PRNT_PTRTYPE "]) from " "%s to %s, saving to %s ...\n",
                params->iFile, docs->nrows, docs->ncols,
                docs->rowptr[docs->nrows], da_getStringKey (fmt_options,
                        params->fmtRead),
                        da_getStringKey (fmt_options, params->fmtWrite), params->oFile);

    // pre-process the data
    preProcessData (params);

    if (params->fmtWrite == DA_FMT_BINCOL && !docs->colptr)
        da_csr_CreateIndex (docs, DA_COL);
    else if (!docs->rowptr)
        da_csr_CreateIndex (docs, DA_ROW);

    da_csr_Write (docs, params->oFile, params->fmtWrite, params->writeVals,
            params->writeNum);
    if (params->verbosity > 0)
        printf ("Done.\n");

    freeParams (&params);
    exit (EXIT_SUCCESS);
}


/**
 * Free memory from the params structure
 */
void freeParams (params_t ** params)
{
    da_csr_FreeAll (&(*params)->docs, &(*params)->docsdb, &(*params)->neighbors, LTERM);
    if ((*params)->fpout)
        da_fclose ((*params)->fpout);
    if ((*params)->neighcache) {
        da_free ((void **) &(*params)->neighcache->sims, &(*params)->neighcache,
                LTERM);
    }
    da_free ((void **) &(*params)->iFile, &(*params)->dFile, &(*params)->oFile,
            &(*params)->vFile, &(*params)->dataset, &(*params)->filename,
            &(*params)->rperm, &(*params)->rsizes, &(*params)->csizes,
            &(*params)->rwgts, &(*params)->cwgts, LTERM);

    da_free ((void **) params, LTERM);
}
