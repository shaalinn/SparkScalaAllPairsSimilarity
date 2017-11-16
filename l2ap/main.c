/*!
 \file  main.c
 \brief This file is the entry point for apss's various components

 \author David C. Anastasiu
 */

#include "includes.h"
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>


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
	sprintf(fname,"/home/shalin/Downloads/l2ap/build/test%ld-%d.clu",now.tv_usec,r);
	

	//This lines for file with header
	sprintf(fname2,"/home/shalin/Downloads/l2ap/build/test%ld-%d(2).clu",now.tv_usec,r);
	
	//Uncomment following two lines to enable doc2mat
	//sprintf(textfile,"/Users/samin/Downloads/l2ap/build/textfile%ld-%d.raw",now.tv_usec,r);
	//sprintf(command,"/Users/samin/Downloads/doc2mat-1.0/./doc2mat -nlskip=1 %s %s >> doc2matop.txt",textfile, fname);
	
	
	
	char str[10000];
	FILE *fp, *fp2,*fp3;
	//printf("hi");
	params = (params_t *)gk_malloc(sizeof(params_t), "main: params");
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
	
	fp3=fopen("/home/shalin/Downloads/l2ap/build/cluHeaders.txt","r");
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
	//pup=fopen("/home/shalin/Downloads/l2ap/build/parth.txt","w+");
	//fputs(threshold,pup);
	//fclose(pup);
	int a=system(command);
	char *args[]={"./apss",threshold,"l2ap",fname2};
	cmdline_parse(params, 4, args);

	if(params->verbosity > 0){
		printf("********************************************************************************\n");
		printf("%s (%d.%d.%d), vInfo: [%s]\n", PROGRAM_NAME, VER_MAJOR, VER_MINOR,
				VER_SUBMINOR, VER_COMMENT);
		//params->mode=params->mode?params->mode:MODE_L2AP;
		//params->sim = DA_SIM_COS;
		printf("mode: %s, ", da_getStringKey(mode_options, params->mode));
		printf("iFile: %s, ", params->iFile? params->iFile : "NULL");
		printf("oFile: %s, ", params->oFile ? params->oFile : "NULL");
		if(params->mode != MODE_TESTEQUAL && params->mode != MODE_IO) {
			printf("t: %.2f", params->simT);
		} else if(params->mode == MODE_TESTEQUAL) {
			printf("fldelta: %g", params->fldelta);
		}
		printf("\n********************************************************************************\n");
		fflush(stdout);
	}

     
	gk_startwctimer(params->timer_global);

	// read input data
	readInputData(params);

	switch(params->mode){

	case MODE_TESTEQUAL:
		da_testMatricesEqual(params);
		break;

	case MODE_INFO:
		da_matrixInfo(params);
		break;

	case MODE_IO:
		da_matrixIo(params);
		break;

	case MODE_IDXJOIN:
			ijFindNeighbors(params);
		break;

	case MODE_AP:
	case MODE_AP2:
			apFindNeighbors(params);
		break;

    case MODE_MMJOIN:
            mmjFindNeighbors(params);
        break;

    case MODE_MKJOIN:
            mkjFindNeighbors(params);
        break;

    case MODE_MKJOIN2:
            mkjFindNeighbors2(params);
        break;

    case MODE_L2AP:
            if(params->sim == DA_SIM_COS){
                l2apFindNeighbors(params);
            } else {
                l2apFindNeighborsTan(params);
            }
        break;

    case MODE_L2AP_T2:
            l2apFindNeighborsTan2(params);
        break;

    case MODE_L2AP_CT:
    case MODE_L2AP_MT:
            l2apFindNeighborsTanM(params);
		break;

	default:
		gk_errexit(SIGERR, "Invalid mode.");
		break;
	}


	// similarity search complete.
	if(params->verbosity > 0){
		printf("Num. total candidates: %zu\n", params->nCandidates);
		printf("Num. dot products: %zu\n", params->nDotProducts);
#ifdef EXTRACOUNTS
		if(params->mode >= MODE_AP){
			printf("Dynamic index size: %zu\n", params->indexSize);
			printf("Nnzs pruned by minsize bound: %zu\n", params->nPruneMinsize);
			if(params->mode >= MODE_L2AP){
				printf("Num. candidates pruned by l2-norm bound (gen): %zu\n", params->nPruneLength);
				printf("Num. candidates pruned by l2-norm bound (ver): %zu\n", params->nPruneLength2);
				printf("Num. candidates pruned by pscore bound: %zu\n", params->nPrunePscore);
			}
			printf("Num. candidates pruned by dotp bound: %zu\n", params->nPruneDotP);
			if(params->mode >= MODE_L2AP){
				printf("Num. candidates pruned by positional dotp bound: %zu\n", params->nPruneDotP2);
				if(params->mode >= MODE_L2AP && params->sim == DA_SIM_TAN){
                    printf("Num. objects for which the Cosine RS bound was effective: %zu\n", params->nCountRS);
                    printf("Num. objects for which the Tanimoto RS bound was effective: %zu\n", params->nCountTanRS);
                    printf("Num. objects ignored due to original vector length: %zu\n", params->nPruneTanLen);
                    printf("Num. candidates pruned in CV due to candidate vector length: %zu\n", params->nPruneTanLenCV);
				}
			}
		} else if(params->mode == MODE_MMJOIN){
			printf("Dynamic index size: %zu\n", params->indexSize);
		} else if(params->mode == MODE_MKJOIN){
            printf("Num. candidates pruned due to original vector length: %zu\n", params->nPruneTanLen);
		}
#endif
		printf("Num. similar pairs: %zu\n", params->nSimPairs);
		if(params->nim){
			printf("Num. times neighbors structure increased in size: " PRNT_PTRTYPE "\n", params->nghinccnt);
		}
		fflush(stdout);

		gk_stopwctimer(params->timer_global);

#ifdef EXTRATIMES
		printf("TIMES:\n");
		da_printTimerLong("\t Add to index: ",
				gk_getwctimer(params->timer_4));
		da_printTimerLong("\t Generate candidates: ",
				gk_getwctimer(params->timer_1));
		da_printTimerLong("\t Process candidates: ",
				gk_getwctimer(params->timer_2));
#endif
		da_printTimerLong("\t Similarity search: ",
				gk_getwctimer(params->timer_3));
		da_printTimerLong("\t Total time: ",
				gk_getwctimer(params->timer_global));

		printf(
				"********************************************************************************\n");
	}

	freeParams(&params);
	//system("mv /Users/samin/Downloads/l2ap/build/test* /Users/samin/Downloads/l2ap/build/untitled_folder/");
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
		gk_startwctimer(params->timer_5); // memory allocation time
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
		params->neighcache = (da_sims_t*)gk_malloc(sizeof(da_sims_t), "l2apFindNeighbors: neighcache");
		params->neighcache->nsims = 0;
		params->neighcache->size = NSIMSBUFFER;
		params->neighcache->sims = da_ssmalloc(NSIMSBUFFER, (da_sim_t){0,0,0.0}, "l2apFindNeighbors: neighcache->sims");

		/* open output file */
		if(params->oFile && params->mode < MODE_INFO
                && gk_strcasecmp(params->oFile, "none") == 0
                && gk_strcasecmp(params->oFile, "null") == 0)
			params->fpout = gk_fopen(params->oFile, "w", "da_csr_Write: fpout");
	}
}



void simSearchFinalize(params_t *params){
	if(params->nim){
		// multiply by transpose to represent similarity with higher rows
		da_addHigherNeighbors(params->neighbors);
		// inverse permute neighbors
		if(params->rperm)
			da_inversePermuteMatrix(&(params->neighbors), params->rperm, params->rperm);
		if(params->verbosity > 0)
			printf("Writing neighborhood matrix to %s.\n", params->oFile);
		if(strncmp(params->oFile, "none", 4) != 0 && strncmp(params->oFile, "null", 4) != 0)
		    da_csr_Write(params->neighbors, params->oFile, params->fmtWrite, 1, 1);
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
		gk_errexit(SIGERR, "Invalid input format.\n");
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
	da_csrCompare(docs, docs2, params->fldelta, 1);
	da_csr_Free(&docs2);
	freeParams(&params);
	exit(EXIT_SUCCESS);
}

/**
 * Display information about a sparse matrix
 */
void da_matrixInfo(params_t *params){
	da_csr_t *docs = params->docs;

	// identify format, read matrix, and print out information about it - nrows, ncols, nnz
	printf("%s: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, " PRNT_PTRTYPE " nnzs, %g density, ",
			params->iFile, docs->nrows, docs->ncols, docs->rowptr[docs->nrows],
			docs->rowptr[docs->nrows] / ((double) docs->nrows * docs->ncols)
	);

	da_csr_CompactColumns(docs);
	printf(PRNT_IDXTYPE " non-empty cols.\n\n", docs->ncols);

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
			gk_errexit(SIGERR, "Invalid output format.\n");
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
		gk_fclose((*params)->fpout);
	if((*params)->neighcache){
		gk_free((void**)&(*params)->neighcache->sims, &(*params)->neighcache, LTERM);
	}
	gk_free((void**)&(*params)->iFile, &(*params)->oFile, &(*params)->dataset,
			&(*params)->filename, &(*params)->rperm, &(*params)->cperm, LTERM);
	gk_free((void**)params, LTERM);
}
