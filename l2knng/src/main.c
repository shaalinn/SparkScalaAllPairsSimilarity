/*!
 \file  main.c
 \brief This file is the entry point for apss's various components

 \author David C. Anastasiu
 */

#include "includes.h"


/***********************************************/
/*! This is the entry point for the program    */
/***********************************************/
int main(int argc, char *argv[]) {
	params_t *params;
	char fname[256],fname2[256],command[256],textfile[256];
	time_t seconds;
	seconds=time(NULL);
	struct timeval now;
	gettimeofday(&now, NULL);

	srand(time(NULL));   // should only be called once
	int r = rand();
	
	//printf("%ld %ld",now.tv_usec,r);
	sprintf(fname,"/home/shalin/Downloads/l2knng/build/test%ld-%d.clu",now.tv_usec,r);
	

	//This lines for file with header
	sprintf(fname2,"/home/shalin/Downloads/l2knng/build/test%ld-%d(2).clu",now.tv_usec,r);
	
	//Uncomment following two lines to enable doc2mat
	//sprintf(textfile,"/Users/samin/Downloads/l2ap/build/textfile%ld-%d.raw",now.tv_usec,r);
	//sprintf(command,"/Users/samin/Downloads/doc2mat-1.0/./doc2mat -nlskip=1 %s %s >> doc2matop.txt",textfile, fname);
	
	
	
	char str[10000];
	FILE *fp, *fp2,*fp3;
	//printf("hi");
	//fp = fopen(textfile, "w+");
	fp=fopen(fname,"w+");
	int rowcount=0;
	while (1) {
		char ch=fgetc(stdin);
	    if (ch==EOF) {	
	      break;
	    }
	    if(ch=='\n'){
	      rowcount++;
	    }
	    fputc(ch,fp);
	 }
	fclose(fp);

	fp2=fopen(fname2,"w+");
	
	fp3=fopen("/home/shalin/Downloads/l2knng/build/cluHeaders.txt","r");
	fprintf(fp2,"%d ",rowcount);
	while(1){
		char ch=fgetc(fp3);
		if(ch=='\n'){
			
			break;
		}
		fputc(ch,fp2);
	}
	fputc('\n',fp2);
	char threshold[20];
	threshold[0]='-';
	threshold[1]='t';
	threshold[2]='=';
	int z=3;
	while(1){
		char ch = fgetc(fp3);
		if(ch==EOF){
			break;
		}
		threshold[z]=ch;
		z++;
	}
	threshold[z]='\0';
	
	fclose(fp3);
	fp=fopen(fname,"r");
	while(1){
		char ch=fgetc(fp);
		if(ch==EOF){
			break;
		}
		fputc(ch,fp2);
	}
	fclose(fp);
	fclose(fp2);

	//FILE *pup;
	//pup=fopen("/home/shalin/Downloads/l2knng/build/parth.txt","w+");
	//fputs(threshold,pup);
	//fclose(pup);
	int a=system(command);
	char *args[]={"./knng",threshold,"l2knn",fname2};
	cmdline_parse(params, 4, args);

	params = (params_t *)da_malloc(sizeof(params_t), "main: params");

	cmdline_parse(params, argc, argv);

	if(params->verbosity > 0){
		printf("********************************************************************************\n");
		printf("%s (%d.%d.%d), vInfo: [%s]\n", PROGRAM_NAME, VER_MAJOR, VER_MINOR,
				VER_SUBMINOR, VER_COMMENT);
		printf("mode: %s, ", da_getStringKey(mode_options, params->mode));
		printf("iFile: %s, ", params->iFile);
		printf("oFile: %s, ", params->oFile ? params->oFile : "NULL");
		if(params->mode == MODE_KIDXJOIN || params->mode == MODE_MSC || params->mode == MODE_BMM || params->mode == MODE_BMMC) {
            printf("k: %d", params->k);
        } else if(params->mode == MODE_KL2AP) {
            printf("k: %d, step: %.3f", params->k, params->step);
        } else if(params->mode == MODE_L2KNN || params->mode == MODE_L2KNN_A) {
            printf("k: %d, alpha: %.4f, enh: %d", params->k, params->alpha, params->nenhance);
        } else if(params->mode == MODE_GF) {
            printf("k: %d, mu: %d", params->k, (int)params->mu);
        } else if(params->mode == MODE_PKIDXJOIN) {
            printf("k: %d, nqrows: %d, ndrows: %d, nthreads: %d", params->k,
                    params->nqrows, params->ndrows, params->nthreads);
        } else if(params->mode == MODE_TESTEQUAL) {
			printf("fldelta: %g", params->fldelta);
		}
		printf("\n********************************************************************************\n");
		fflush(stdout);
	}

	da_startwctimer(params->timer_global);

	// read input data
	readInputData(params);

	switch(params->mode){

	case MODE_KL2AP:
        kl2apFindNeighbors(params);
		break;

    case MODE_KIDXJOIN:
        kijFindNeighbors(params);
        break;

    case MODE_PKIDXJOIN:
        pkijFindNeighbors(params);
        break;

    case MODE_L2KNN:
    case MODE_L2KNN_A:
        l2knnFindNeighbors(params);
        break;

    case MODE_MSC:
        msFindNeighbors(params);
        break;

    case MODE_BMM:
        bmmFindNeighbors(params);
        break;

    case MODE_BMMC:
        bmmcFindNeighbors(params);
        break;

    case MODE_GF:
        gfFindNeighbors(params);
        break;



    case MODE_TESTEQUAL:
        da_testMatricesEqual(params);
        break;

    case MODE_INFO:
        da_matrixInfo(params);
        break;

    case MODE_IO:
        da_matrixIo(params);
        break;

    case MODE_RECALL:
        da_testRecall(params);
        break;

	default:
		da_errexit("Invalid mode.");
		break;
	}


	// similarity search complete.
	if(params->verbosity > 0){
		printf("Num. total candidates: %zu\n", params->nCandidates);
		printf("Num. dot products: %zu\n", params->nDotProducts);
#ifdef EXTRACOUNTS
        if(params->mode == MODE_L2AP){
            printf("Dynamic index size: %zu\n", params->indexSize);
            printf("Nnzs pruned by minsize bound: %zu\n", params->nPruneMinsize);
            printf("Num. candidates pruned by pscore bound: %zu\n", params->nPrunePscore);
            printf("Num. candidates pruned by dotp bound: %zu\n", params->nPruneDotP);
            printf("Num. candidates pruned by positional dotp bound: %zu\n", params->nPruneDotP2);
        }
        if(params->mode == MODE_L2KNN || params->mode == MODE_L2KNN_A){
            printf("Num. initial construction dot products: %zu\n", params->nDotProducts1);
            printf("Num. enhancement dot products: %zu\n", params->nDotProducts2);
            printf("Num. candidates pruned by l2-norm bound: %zu\n", params->nPruneLength);
            printf("Num. candidates pruned by pscore bound: %zu\n", params->nPrunePscore);
            printf("Num. candidates pruned during verification by l2-norm bound: %zu\n", params->nPruneLength2);
            printf("Dynamic index size: %zu\n", params->indexSize);
        }

#endif
		printf("Num. similar pairs: %zu\n", params->nSimPairs);
		if(params->nim && params->mode == MODE_L2AP)
		    printf("Num. times neighbors structure increased in size: " PRNT_PTRTYPE "\n", params->nghinccnt);
		fflush(stdout);

		da_stopwctimer(params->timer_global);

#ifdef EXTRATIMES
		printf("TIMES:\n");
		da_printTimerLong("\t Add to index: ",
				da_getwctimer(params->timer_4));
		da_printTimerLong("\t Generate candidates: ",
				da_getwctimer(params->timer_1));
		da_printTimerLong("\t Process candidates: ",
				da_getwctimer(params->timer_2));
#endif
		da_printTimerLong("\t Similarity search: ",
				da_getwctimer(params->timer_3));
		da_printTimerLong("\t Total time: ",
				da_getwctimer(params->timer_global));

		printf(
				"********************************************************************************\n");
	}

	freeParams(&params);
	exit(EXIT_SUCCESS);
}

/**
 * Pre-process dataset by pruning, scaling, filtering, shuffling, etc.
 */
void preProcessData(params_t *params){
	da_csr_t *docs = params->docs;
	da_csr_t *tmp = NULL;

	/* prune cols with less than pcminlen or more than pcmaxlen values */
	if(params->prunecols){
		if(params->pcminlen == -1)
			params->pcminlen = 0;
		if(params->pcmaxlen == -1)
			params->pcmaxlen = docs->nrows;
		printf("pc: %d, pcmin: %d, pcmax: %d\n", params->prunecols, params->pcminlen, params->pcmaxlen);
		tmp = da_csr_Prune(docs, DA_COL, params->pcminlen, params->pcmaxlen);
		da_csr_Transfer(tmp, docs);
		da_csr_Free(&tmp);
	}

	/* prune rows with less than prminlen or more than prmaxlen values */
	if(params->prunerows){
		if(params->prminlen == -1)
			params->prminlen = 0;
		if(params->prmaxlen == -1)
			params->prmaxlen = docs->ncols;
		printf("pr: %d, prmin: %d, prmax: %d\n", params->prunerows, params->prminlen, params->prmaxlen);
		tmp = da_csr_Prune(docs, DA_ROW, params->prminlen, params->prmaxlen);
		da_csr_Transfer(tmp, docs);
		da_csr_Free(&tmp);
	}


    // compact the column space
    if(params->compactcols)
        da_csr_CompactColumns(docs);

    // sort the column space
    if(params->sortcols)
        da_csr_SortIndices(docs, DA_ROW);

	// compact the row space
	if(params->compactrows)
		da_csr_CompactRows(docs);

	/* scale term values */
	if(params->scale){
		if(params->verbosity > 0)
			printf("   Scaling input matrix.\n");
		da_csr_Scale(docs, DA_SCALE_IDF);
	}

	/* normalize docs rows */
	if(params->norm > 0)
		da_csr_Normalize(docs, DA_ROWCOL, params->norm);

}

void simSearchSetup(params_t *params){
	da_csr_t *docs = params->docs;
	da_csr_t *neighbors = NULL;
	size_t nghnnz;

	// pre-process the data
	preProcessData(params);

	// compact the column space
	da_csr_CompactColumns(docs);

	if(params->verbosity > 0)
		printf("Docs matrix: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, "
			PRNT_PTRTYPE " nnz\n", docs->nrows, docs->ncols, docs->rowptr[docs->nrows]);

	if(params->nim){
		// allocate memory for tracking neighbors - O <-- {}
		da_startwctimer(params->timer_5); // memory allocation time
		neighbors = da_csr_Create();
		neighbors->nrows = neighbors->ncols = docs->nrows;
		nghnnz = params->nghnnz = NINITNEIGHBORS * docs->nrows;
		neighbors->rowptr = da_pmalloc(docs->nrows + 1, "apFindNeighbors: neighbors->rowptr");
		neighbors->rowind = da_imalloc(nghnnz, "apFindNeighbors: neighbors->rowind");
		neighbors->rowval = da_vmalloc(nghnnz, "apFindNeighbors: neighbors->rowval");
		neighbors->rowptr[0] = 0;
		params->neighbors = neighbors;

	} else {
		// set up neighbor output cache
		params->neighcache = (da_sims_t*)da_malloc(sizeof(da_sims_t), "l2apFindNeighbors: neighcache");
		params->neighcache->nsims = 0;
		params->neighcache->size = NSIMSBUFFER;
		params->neighcache->sims = da_ssmalloc(NSIMSBUFFER, (da_sim_t){0,0,0.0}, "l2apFindNeighbors: neighcache->sims");

		/* open output file */
		if(params->oFile && params->mode < MODE_INFO)
			params->fpout = da_fopen(params->oFile, "w", "da_csr_Write: fpout");
	}
}



void simSearchFinalize(params_t *params){
	if(params->nim){
		// multiply by transpose to represent similarity with higher rows
		da_addHigherNeighbors(params->neighbors);
		// inverse permute neighbors
		if(params->rperm)
			da_inversePermuteMatrix(&(params->neighbors), params->rperm, params->rperm);
		if(params->fpout) {
            if(params->verbosity > 0)
                printf("Writing neighborhood matrix to %s.\n", params->oFile);
            da_csr_Write(params->neighbors, params->oFile, params->fmtWrite, 1, 1);
		}
	} else {
		// write out neighbors buffer
		da_storeSim(params, -1, -1, 0);
	}
}

/**
 * Read input data
 */
void readInputData(params_t *params){
	da_csr_t *docs;
	params->fmtRead = da_getFileFormat(params->iFile, params->fmtRead);
	if(params->fmtRead < 1)
		da_errexit("Invalid input format.\n");
	docs = da_csr_Read(params->iFile, params->fmtRead, params->readVals, params->readNum);
	assert(docs->rowptr || docs->colptr);
	/* ensure row based structure for sim search */
	if(params->mode < 50 && !docs->rowptr){
		da_csr_CreateIndex(docs, DA_ROW);
		da_csr_FreeBase(docs, DA_COL);
	}
	params->docs = docs;
}

/**
 * Test two matrices are equal. Values tested up to params->fldelta precision.
 */
void da_testMatricesEqual(params_t *params){
	da_csr_t *docs=NULL, *docs2=NULL;

	docs = params->docs;

	// test equality of two sparse matrices & print out differences
	docs2 = da_csr_Read(params->oFile, da_getFileFormat(params->oFile, params->fmtWrite),
			params->readVals, params->readNum);
	printf("Comparing %s (A[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE "]) and "
			"%s (B[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE "]).\n\n",
			params->iFile, docs->nrows, docs->ncols, docs->rowptr[docs->nrows],
			params->oFile, docs2->nrows, docs2->ncols, docs2->rowptr[docs2->nrows]);
	da_csrCompare(docs, docs2, params->fldelta, params->cmpind, params->cmpval);
	da_csr_Free(&docs2);
	freeParams(&params);
	exit(EXIT_SUCCESS);
}



/**
 * Test two matrices are equal. Values tested up to params->fldelta precision.
 */
void da_testRecall(params_t *params){
    da_csr_t *docs=NULL, *docs2=NULL;

    docs = params->docs;

    // test equality of two sparse matrices & print out differences
    docs2 = da_csr_Read(params->oFile, da_getFileFormat(params->oFile, params->fmtWrite),
            params->readVals, params->readNum);
    printf("Usage: knng recall <test_results> <true_results>\n"
            "Use -verb 3 for additional information. Neighbors will be marked with:\n"
            "\t* neighbors we missed with same value as the min values\n"
            "\t+ neighbors we reported that are not in the true neighborhood\n"
            "\t- neighbors we missed and did not report\n\n");
    printf("Comparing %s (A[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE "]) and "
            "%s (B[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE "]).\n\n",
            params->iFile, docs->nrows, docs->ncols, docs->rowptr[docs->nrows],
            params->oFile, docs2->nrows, docs2->ncols, docs2->rowptr[docs2->nrows]);
    verify_knng_results3(params->docs, docs2, params->k, params->verbosity);
    da_csr_Free(&docs2);
    freeParams(&params);
    exit(EXIT_SUCCESS);
}

/**
 * Display information about a sparse matrix
 */
void da_matrixInfo(params_t *params){
	da_csr_t *docs = params->docs;
	size_t i, j, l, x1, x2, min, max;

	// identify format, read matrix, and print out information about it - nrows, ncols, nnz
	printf("%s: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, " PRNT_PTRTYPE " nnzs, %g density, ",
			params->iFile, docs->nrows, docs->ncols, docs->rowptr[docs->nrows],
			docs->rowptr[docs->nrows] / ((double) docs->nrows * docs->ncols)
	);

	da_csr_CompactColumns(docs);
	printf(PRNT_IDXTYPE " non-empty cols.\n", docs->ncols);

	if(params->stats){
	    if(!docs->colptr)
	        da_csr_CreateIndex(docs, DA_COL);
	    for(x1=0, x2=0, min=INT_MAX, max=0, i=0; i < docs->nrows; ++i){
	        l = docs->rowptr[i+1] - docs->rowptr[i];
	        x1 += l;
	        x2 += l*l;
	        if(l > max)
	            max = l;
	        if(l < min)
	            min = l;
	    }
	    printf("Row nnz stats: min %zu, max %zu mean %.2f, stdev %.2f.\n", min, max,
	            (double)x1/(double)docs->nrows,
	            sqrt( (double)x2*(double)docs->nrows - x1*x1 ) / (double)docs->nrows );
	    for(x1=0, x2=0, i=0; i < docs->ncols; ++i){
            l = docs->colptr[i+1] - docs->colptr[i];
            x1 += l;
            x2 += l*l;
            if(l > max)
                max = l;
            if(l < min)
                min = l;
        }
        printf("Col nnz stats: min %zu, max %zu mean %.2f, stdev %.2f.\n", min, max,
                (double)x1/(double)docs->ncols,
                sqrt( (double)x2*(double)docs->ncols - x1*x1 ) / (double)docs->ncols );
	}



	printf("\n");

	freeParams(&params);
	exit(EXIT_SUCCESS);
}


/**
 * Transform input matrix to some other format
 * Data pre-processing invoked before transforming.
 */
void da_matrixIo(params_t *params){
	da_csr_t *docs = params->docs;

	if(params->oFile){
		params->fmtWrite = da_getFileFormat(params->oFile, params->fmtWrite);
		if(params->fmtWrite < 1)
			da_errexit("Invalid output format.\n");
	}

	if(params->fmtRead == DA_FMT_BINAPB && params->fmtWrite != DA_FMT_BINAPB){
		// transition from binary to weighted format
		if(!docs->rowval)
			docs->rowval = da_vsmalloc(docs->rowptr[docs->nrows], 1.0, "io: docs->rowptr");
	}

	if(params->verbosity > 0)
		printf("Transforming %s (A[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE "]) from "
			"%s to %s, saving to %s ...\n",
			params->iFile, docs->nrows, docs->ncols, docs->rowptr[docs->nrows],
			da_getStringKey(fmt_options, params->fmtRead),
			da_getStringKey(fmt_options, params->fmtWrite), params->oFile);

	// pre-process the data
	preProcessData(params);

	if(params->fmtWrite == DA_FMT_BINCOL && ! docs->colptr)
		da_csr_CreateIndex(docs, DA_COL);
	else if(!docs->rowptr)
		da_csr_CreateIndex(docs, DA_ROW);

	da_csr_Write(docs, params->oFile, params->fmtWrite, params->writeVals, params->writeNum);
	if(params->verbosity > 0)
		printf("Done.\n");

	freeParams(&params);
	exit(EXIT_SUCCESS);
}


/**
 * Free memory from the params structure
 */
void freeParams(params_t** params){
	da_csr_FreeAll(&(*params)->docs, &(*params)->neighbors, LTERM);
	if((*params)->fpout)
		da_fclose((*params)->fpout);
	if((*params)->neighcache){
		da_free((void**)&(*params)->neighcache->sims, &(*params)->neighcache, LTERM);
	}
	da_free((void**)&(*params)->iFile, &(*params)->oFile, &(*params)->dataset, &(*params)->vFile,
	        &(*params)->filename, &(*params)->crfname, &(*params)->crfname2, &(*params)->querydocs, &(*params)->rperm,
			&(*params)->rsizes, &(*params)->csizes, &(*params)->rwgts, &(*params)->cwgts,
			LTERM);

	da_free((void**)params, LTERM);
}
