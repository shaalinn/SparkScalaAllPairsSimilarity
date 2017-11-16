/*!
 \file  util.c
 \brief This file contains utilities for the application

 \author David C. Anastasiu
 */

#include "includes.h"

/** forward declarations **/
int da_crRead2(params_t *params, da_iapq_t **knng, val_t *simT, size_t *rid,
        FILE* fpin, char* f_str, ...);


/*************************************************************************/
/*! This function prints an error message and raises a signum signal
 */
/*************************************************************************/
void da_errexit(const char* const f_str,...)
{
    va_list argp;

    va_start(argp, f_str);
    vfprintf(stderr, f_str, argp);
    va_end(argp);

    fprintf(stderr,"\n");
    fflush(stderr);

    raise(SIGTERM);
}




/*************************************************************************
* This function returns the CPU seconds
**************************************************************************/
double da_WClockSeconds(void)
{
#ifdef __GNUC__
    struct timeval ctime;

    gettimeofday(&ctime, NULL);

    return (double)ctime.tv_sec + (double).000001*ctime.tv_usec;
#else
    return (double)time(NULL);
#endif
}


/*************************************************************************
 * This function returns the CPU seconds
 **************************************************************************/
double da_CPUSeconds(void)
{
    //#ifdef __OPENMP__
#ifdef __OPENMPXXXX__
    return omp_get_wtime();
#else
#if defined(WIN32) || defined(__MINGW32__)
    return((double) clock()/CLOCKS_PER_SEC);
#else
    struct rusage r;

    getrusage(RUSAGE_SELF, &r);
    return ((r.ru_utime.tv_sec + r.ru_stime.tv_sec) + 1.0e-6*(r.ru_utime.tv_usec + r.ru_stime.tv_usec));
#endif
#endif
}

/**
 * Write out a checkpoint restart for the knng solution in kIdxJoin and kl2ap modes
 */
void da_crWrite(params_t *params, da_iapq_t **knng, val_t simT, size_t rid)
{
    size_t i, j, type;
    idx_t nrows;
    ptr_t nnz;
    FILE* fpout;

    nrows = params->docs->nrows;

    if(!params->crfname)
        da_errexit("da_crWrite: Checkpoint restart filename required for checkpoint restart!");

    /* move backup solution if it exists */
    if(da_fexists(params->crfname)){
        unlink(params->crfname2);
        if(rename(params->crfname, params->crfname2) != 0)
            da_errexit("Could not rename existing cr file!");
    }

    /* stop active timers */
    da_stopwctimer(params->timer_3);
    da_stopwctimer(params->timer_global);

    fpout = da_fopen(params->crfname, "w", "da_vWriteMatrix: fpout");

    /* store params book-keeping data */
    fwrite(&(params->nCandidates), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nDotProducts), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nDotProducts1), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nDotProducts2), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nPruneDotP), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nPruneDotP2), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nPruneLength), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nPruneLength2), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nPruneMinsize), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nPrunePscore), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->nSimPairs), sizeof(ssize_t), 1, fpout);
    fwrite(&(params->indexSize), sizeof(ssize_t), 1, fpout);

    /* store timers */
    fwrite(&(params->timer_global), sizeof(double), 1, fpout);
    fwrite(&(params->timer_1), sizeof(double), 1, fpout);
    fwrite(&(params->timer_2), sizeof(double), 1, fpout);
    fwrite(&(params->timer_3), sizeof(double), 1, fpout);
    fwrite(&(params->timer_4), sizeof(double), 1, fpout);
    fwrite(&(params->timer_5), sizeof(double), 1, fpout);
    fwrite(&(params->timer_6), sizeof(double), 1, fpout);
    fwrite(&(params->timer_7), sizeof(double), 1, fpout);
    fwrite(&(params->timer_8), sizeof(double), 1, fpout);
    fwrite(&(params->timer_9), sizeof(double), 1, fpout);

    /* store simT and rid */
    fwrite(&simT, sizeof(val_t), 1, fpout);
    fwrite(&rid, sizeof(size_t), 1, fpout);

    if(knng){
        /* store the current version of the knng */
        type = 1;
        fwrite(&type, sizeof(size_t), 1, fpout);
        fwrite(&nrows, sizeof(idx_t), 1, fpout); /* number of knn lists */
        for(i=0; i < nrows; ++i){
            fwrite(&(knng[i]->nnodes), sizeof(size_t), 1, fpout); /* length of the list */
            assert(knng[i]->nnodes <= params->k);
            if(knng[i]->nnodes > 0)
                fwrite(knng[i]->heap, sizeof(da_iakv_t), knng[i]->nnodes, fpout); /* heap list */
        }
    } else if(params->neighbors) {
        /* store partial neighbors CSR */
        type = 2;
        fwrite(&type, sizeof(size_t), 1, fpout);
        nnz = params->neighbors->rowptr[rid+1];
        fwrite(&nnz, sizeof(ptr_t), 1, fpout);
        fwrite(params->neighbors->rowptr, sizeof(ptr_t), rid+1, fpout);
        fwrite(params->neighbors->rowind, sizeof(idx_t), nnz, fpout);
        fwrite(params->neighbors->rowval, sizeof(val_t), nnz, fpout);
    } else
        da_errexit("da_crWrite: Both knng and neighbors missing!");

    da_fclose(fpout);

    da_startwctimer(params->timer_3);
    da_startwctimer(params->timer_global);

}


/**
 * Read in a checkpoint restart for the knng solution in kIdxJoin and kl2ap modes
 */
int da_crRead(params_t *params, da_iapq_t **knng, val_t *simT, size_t *rid)
{
    size_t i, j, len, id, type;
    idx_t nrows;
    val_t st;
    ptr_t nnz;
    FILE* fpin;

    if(!params->crfname)
        da_errexit("da_crRead: Checkpoint restart filename required for checkpoint restart!");


    da_stopwctimer(params->timer_3);
    da_stopwctimer(params->timer_global);

    fpin = da_fopen(params->crfname, "rb", "da_csr_Read: fpin");

    /* read in the params data */
    if (fread(&(params->nCandidates), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nCandidates from file %s!\n", params->crfname);
    if (fread(&(params->nDotProducts), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nDotProducts from file %s!\n", params->crfname);
    if (fread(&(params->nDotProducts1), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nDotProducts1 from file %s!\n", params->crfname);
    if (fread(&(params->nDotProducts2), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nDotProducts2 from file %s!\n", params->crfname);
    if (fread(&(params->nPruneDotP), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nPruneDotP from file %s!\n", params->crfname);
    if (fread(&(params->nPruneDotP2), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nPruneDotP2 from file %s!\n", params->crfname);
    if (fread(&(params->nPruneLength), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nPruneLength from file %s!\n", params->crfname);
    if (fread(&(params->nPruneLength2), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nPruneLength2 from file %s!\n", params->crfname);
    if (fread(&(params->nPruneMinsize), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nPruneMinsize from file %s!\n", params->crfname);
    if (fread(&(params->nPrunePscore), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nPrunePscore from file %s!\n", params->crfname);
    if (fread(&(params->nSimPairs), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read nSimPairs from file %s!\n", params->crfname);
    if (fread(&(params->indexSize), sizeof(ssize_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read indexSize from file %s!\n", params->crfname);

    /* read in the timers */
    if (fread(&(params->timer_global), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_global from file %s!\n", params->crfname);
    if (fread(&(params->timer_1), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_1 from file %s!\n", params->crfname);
    if (fread(&(params->timer_2), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_2 from file %s!\n", params->crfname);
    if (fread(&(params->timer_3), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_3 from file %s!\n", params->crfname);
    if (fread(&(params->timer_4), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_4 from file %s!\n", params->crfname);
    if (fread(&(params->timer_5), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_5 from file %s!\n", params->crfname);
    if (fread(&(params->timer_6), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_6 from file %s!\n", params->crfname);
    if (fread(&(params->timer_7), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_7 from file %s!\n", params->crfname);
    if (fread(&(params->timer_8), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_8 from file %s!\n", params->crfname);
    if (fread(&(params->timer_9), sizeof(double), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read timer_9 from file %s!\n", params->crfname);

    /* read in simT and rid */
    if (fread(&st, sizeof(val_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read simT from file %s!\n", params->crfname);
    if (fread(&id, sizeof(size_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read rid from file %s!\n", params->crfname);
    if (fread(&type, sizeof(size_t), 1, fpin) != 1)
        return da_crRead2(params, knng, simT, rid, fpin,
            "Failed to read type from file %s!\n", params->crfname);


    /* read in the knng */
    if(knng && type == 1){
        if (fread(&nrows, sizeof(idx_t), 1, fpin) != 1)
            return da_crRead2(params, knng, simT, rid, fpin,
                "Failed to read nrows from file %s!\n", params->crfname);
        if(params->docs->nrows != nrows)
            return da_crRead2(params, knng, simT, rid, fpin,
                "Invalid nrows read from file (%d) does not match input nrows (%d)!",
                nrows, params->docs->nrows);
        for(i=0; i < nrows; ++i){
            if (fread(&len, sizeof(size_t), 1, fpin) != 1)
                return da_crRead2(params, knng, simT, rid, fpin,
                    "Failed to read the knng heap length for heap %d from file %s!\n", i, params->crfname);
            knng[i]->nnodes = len;
            if(len > 0){
                ASSERT(len <= knng[i]->maxsize)
                if (fread(knng[i]->heap, sizeof(da_iakv_t), len, fpin) != len)
                    return da_crRead2(params, knng, simT, rid, fpin,
                        "Failed to read the heap for knng heap %d from file %s!\n", i, params->crfname);
            }
        }
    } else if(params->neighbors && type == 2){
        /* read in partial neighbors CSR */
        if (fread(&(nnz), sizeof(ptr_t), 1, fpin) != 1)
            da_errexit( "Failed to read nnz from file %s!\n", params->crfname);
        ASSERT(nnz <= params->nghnnz)
        ASSERT(id <= params->neighbors->nrows)
        if (fread(params->neighbors->rowptr, sizeof(ptr_t), id+1, fpin) != id+1)
            return da_crRead2(params, knng, simT, rid, fpin,
                "Failed to read neighbors rowptr from file %s!\n", params->crfname);
        if (fread(params->neighbors->rowind, sizeof(idx_t), nnz, fpin) != nnz)
            return da_crRead2(params, knng, simT, rid, fpin,
                "Failed to read neighbors rowind from file %s!\n", params->crfname);
        if (fread(params->neighbors->rowval, sizeof(val_t), nnz, fpin) != nnz)
            return da_crRead2(params, knng, simT, rid, fpin,
                "Failed to read neighbors rowval from file %s!\n", params->crfname);
    } else
        return da_crRead2(params, knng, simT, rid, fpin,
            "da_crRead: Both knng and neighbors missing or wrong type: %zu!", type);

    if(simT != NULL)
        *simT = st;
    if(rid != NULL)
        *rid = id;

    da_fclose(fpin);

    printf("Checkpoint restart from %s.\n", params->crfname);

    da_startwctimer(params->timer_3);
    da_startwctimer(params->timer_global);

    return 0;
}

/**
 * In case of checkpoint restart error, try to read from backup
 */
int da_crRead2(params_t *params, da_iapq_t **knng, val_t *simT, size_t *rid,
        FILE* fpin, char* f_str, ...)
{
    if(params->crfname2 && da_fexists(params->crfname2)){
        char *fname;
        int ret;
        da_fclose(fpin);
        fname = params->crfname;
        params->crfname = params->crfname2;
        params->crfname2 = NULL;
        ret = da_crRead(params, knng, simT, rid);
        params->crfname2 = params->crfname;
        params->crfname = fname;

        return ret;
    }
    va_list argp;

    va_start(argp, f_str);
    vfprintf(stderr, f_str, argp);
    va_end(argp);

    fprintf(stderr,"\n");
    fflush(stderr);

    raise(SIGTERM);
    return 1;
}

/*************************************************************************/
/*! If format not specifically given (> 0), check if a text (non-binary) text file
 *  containing a csr is in CLUTO or CSR format.
    \param file is the matrix file to be checked.
    \return the CSR format: DA_FMT_CLUTO or DA_FMT_CSR
 */
/*************************************************************************/
char da_getFileFormat(char *file, const char format)
{
	if(format > 0) return format;
	size_t nnz;
	char fmt;
	char *ext, *p;

	ext = strrchr(file, '.');
	if(ext){
		ext++;
		//make lowercase
		for (p=ext ; *p; ++p) *p = tolower(*p);
		if ((fmt = da_getStringID(fmt_options, ext)) > -1)
			return fmt;
	} else if(da_fexists(file)){ // assume some sort of CSR. Can we guess?
		da_getfilestats(file, NULL, &nnz, NULL, NULL);
		return (nnz%2 == 1) ? DA_FMT_CLUTO : DA_FMT_CSR;
	}
	return -1;
}

/**
 * Retrieve the dataset from the input filename
 */
char* da_getDataset(params_t *params){
	if(params->dataset)
		return params->dataset;
	int len;
	char *end, *start = strrchr(params->iFile, '/');
	if(start)
		start++;
	else
		start = params->iFile;
	end = strstr(start, ".");
	if(!end)
		end = params->iFile + strlen(params->iFile);
	len = end - start;
	params->dataset = da_cmalloc(len+1, "da_getDataset: dataset");
	strncpy(params->dataset, start, len);
	params->dataset[len] = '\0';
	return params->dataset;
}


void da_vWriteVector(char* filename, val_t* vec, idx_t size, char *separator){
	idx_t i;
	FILE* fpout;

	if (filename)
		fpout = da_fopen(filename, "w", "da_vWriteVector: fpout");
	else
		fpout = stdout;

	for(i=0; i < size; i++)
		fprintf(fpout, PRNT_VALTYPE "%s", vec[i], separator);

	if (filename)
		da_fclose(fpout);
}

void da_vWriteMatrix(char* filename, val_t** mat, idx_t nrows, idx_t ncols){
	idx_t i, j;
	FILE* fpout;

	if (filename)
		fpout = da_fopen(filename, "w", "da_vWriteMatrix: fpout");
	else
		fpout = stdout;

	for(i=0; i < nrows; i++){
		for(j=0; j < ncols-1; j++)
			fprintf(fpout, PRNT_VALTYPE "\t", mat[i][j]);
		if(j < ncols)
			fprintf(fpout, PRNT_VALTYPE, mat[i][j]);
		fprintf(fpout, "\n");
	}

	if (filename)
		da_fclose(fpout);
}

void da_dWriteMatrix(char* filename, double** mat, idx_t nrows, idx_t ncols){
	idx_t i, j;
	FILE* fpout;

	if (filename)
		fpout = da_fopen(filename, "w", "da_vWriteMatrix: fpout");
	else
		fpout = stdout;

	for(i=0; i < nrows; i++){
		for(j=0; j < ncols-1; j++)
			fprintf(fpout, "%f\t", mat[i][j]);
		if(j < ncols)
			fprintf(fpout, "%f", mat[i][j]);
		fprintf(fpout, "\n");
	}

	if (filename)
		da_fclose(fpout);
}

void da_iWriteVector(char* filename, idx_t* vec, idx_t size, char *separator){
	idx_t i;
	FILE* fpout;

	if (filename)
		fpout = da_fopen(filename, "w", "da_iWriteVector: fpout");
	else
		fpout = stdout;

	for(i=0; i < size; i++)
		fprintf(fpout, PRNT_IDXTYPE "%s", vec[i], separator);

	if (filename)
		da_fclose(fpout);
}


void da_pWriteVector(char* filename, idx_t* vec, idx_t size, char *separator){
	idx_t i;
	FILE* fpout;

	if (filename)
		fpout = da_fopen(filename, "w", "da_iWriteVector: fpout");
	else
		fpout = stdout;

	for(i=0; i < size; i++)
		fprintf(fpout, PRNT_IDXTYPE "%s", vec[i], separator);

	if (filename)
		da_fclose(fpout);
}

/*************************************************************************/
/*! Stores similarities in the neighbors matrix or writes them out to a file
 * in ijv format. Uses cache to buffer writes. Empties buffer when i is negative.
 *
 *	\param fpout file to write output to
 *	\param cache structure to hold similarity values temporarily
 *	\param zidx Starting col index identifier (1 if first col is 1)
 *	\param i the query neighbor
 *	\param j the other of the two neighbors
 *	\param v similarity value
 */
/*************************************************************************/
void da_storeSim(params_t *params, idx_t i, idx_t j, val_t v){

	FILE *fpout;
	da_sims_t *cache;
	idx_t zidx;
	int32_t k;
	da_sim_t *sims;
	ssize_t nneighbs;

	if(params->nim){  /* we write neighbors in memory for sake of experiments, even if not saving to stdout or disk */
		nneighbs = params->neighbors->rowptr[i+1];
		params->neighbors->rowind[nneighbs] = j;
		params->neighbors->rowval[nneighbs] = v;
		params->neighbors->rowptr[i+1]++;
	} else if(params->fpout) {

		fpout = params->fpout;
		cache = params->neighcache;
		zidx  = params->writeNum;

		if(i >= 0){
			if(params->rperm){
				i = params->rperm[i];
				j = params->rperm[j];
			}
			cache->sims[cache->nsims].i = i;
			cache->sims[cache->nsims].j = j;
			cache->sims[cache->nsims++].val = v;
		}

		if(i < 0 || cache->nsims == cache->size){
			sims = cache->sims;
			for(k=0; k < cache->nsims; k++){
				fprintf(fpout, "%d %d %.8f\n%d %d %.8f\n", (int32_t) sims[k].i + zidx,
						(int32_t) sims[k].j + zidx, (float) sims[k].val, (int32_t) sims[k].j + zidx,
						(int32_t) sims[k].i + zidx, (float) sims[k].val);
			}
			cache->nsims = 0;
		}

	}
}

/**
 * Store sim value in KNNG
 */
void da_storeSimKnn(params_t *params, idx_t i, idx_t j, val_t val)
{
    idx_t rid, cid;
    da_iapq_t **knng = params->knng;

    rid = params->rperm[i];
    cid = params->rperm[j];

    if(knng[rid]->nnodes < params->k || val > da_iapqSeeTopVal(knng[rid]))
        da_iapqInsertHeap(knng[rid], cid, val);
    if(knng[cid]->nnodes < params->k || val > da_iapqSeeTopVal(knng[cid]))
        da_iapqInsertHeap(knng[cid], rid, val);

}


char da_isFmtBinary(char fmt){
	switch(fmt){
	case DA_FMT_BIJV:
	case DA_FMT_BINAP:
	case DA_FMT_BINAPB:
	case DA_FMT_BINCOL:
	case DA_FMT_BINROW:
		return 1;
	default:
		return 0;
	}
}


/*************************************************************************/
/*! Prints a timer in human readable format. Timers should be started and
 * stopped with da_startwctimer(t) and da_stopwctimer(t) respectively.
 *
 *	\param name is the timer name
 *	\param time is the time in seconds
 */
/*************************************************************************/
void da_printTimer(char *name, double time){
	if(time < 60){
		printf("%s %.2fs\n", name, time);
	} else if(time < 3600){
		int32_t min = time/60;
		time -= min*60;
		printf("%s %dm %.2fs\n", name, min, time);
	} else {
		int32_t hours = time/3600;
		time -= hours*3600;
		int32_t min = time/60;
		time -= min*60;
		printf("%s %dh %dm %.2fs\n", name, hours, min, time);
	}
}

/************************************************************************/
/*! \brief Prints a timer in long format, including human readable and
 * computer readable time. Timers should be started and stopped with
 * da_startwctimer(t) and da_stopwctimer(t) respectively.

 	\param name is the Timer name, a string to be printed before the timer.
 	\param time is the time that should be reported.
*/
/*************************************************************************/
void da_printTimerLong(char *name, double time){
	double t = time;
	if(time < 60){
		printf("%s %.2f (%.2fs)\n", name, t, time);
	} else if(time < 3600){
		int32_t min = time/60;
		time -= min*60;
		printf("%s %.2f (%dm %.2fs)\n", name, t, min, time);
	} else {
		int32_t hours = time/3600;
		time -= hours*3600;
		int32_t min = time/60;
		time -= min*60;
		printf("%s %.2f (%dh %dm %.2fs)\n", name, t, hours, min, time);
	}
}


/**
 * http://www.codemaestro.com/reviews/9
 */
float da_fsqrt(float number) {
    long i;
    float x, y;
    const float f = 1.5F;

    x = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;
    i  = 0x5f3759df - ( i >> 1 );
    y  = * ( float * ) &i;
    y  = y * ( f - ( x * y * y ) );
    y  = y * ( f - ( x * y * y ) );
    return number * y;
}

float da_finvSqrt(float x)
{
  float xhalf = 0.5f*x;
  int i = *(int*)&x; // get bits for floating value
  i = 0x5f375a86- (i>>1); // gives initial guess y0
  x = *(float*)&i; // convert bits back to float
  x = x*(1.5f-xhalf*x*x); // Newton step, repeating increases accuracy
  return x;
}



/**
 * Add the transpose of a matrix to its row structure, ignoring self-similarity
 * NOTE1: Assumes rowind has sorted columns
 * NOTE2: Method allocates nnz*2 memory for nval & nrowval.
 *        New matrix may have less nnzs. One may wish to check and
 *        realloc these arrays to free unused space. We leave this
 *        as an optional external procedure.
 */
void da_addHigherNeighbors(da_csr_t* mat){
	ssize_t i, j, k, l, nnz, nrows;
	ptr_t *nrowptr = NULL, *rowptr = NULL, *colptr = NULL;
	idx_t *rowind = NULL, *colind = NULL, *nrowind = NULL, *itmp = NULL;
	val_t *nrowval = NULL, *rowval = NULL, *colval = NULL, *ftmp = NULL;

	ASSERT(mat->nrows == mat->ncols);

	if(mat->nrows < 1)
		return;

	nnz = mat->rowptr[mat->nrows];
	nrows = mat->nrows;
	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;
	nrowptr = da_pmalloc(nrows+1, "da_addHigherNeighbors: nrowptr");
	nrowind = da_imalloc(nnz*2, "da_addHigherNeighbors: nrowind");
	nrowval = da_vmalloc(nnz*2, "da_addHigherNeighbors: nrowval");
	nrowptr[0] = 0;

	// create the transpose of the matrix in the col strucure of mat
	da_csr_CreateIndex(mat, DA_COL);
	da_csr_SortIndices(mat, DA_COL);
	// by definition, the col structure of the mat is equivalent to the transpose of the
	// matrix if thinking of it as row structure instead of col structure
	colptr = mat->colptr;
	colind = mat->colind;
	colval = mat->colval;

	// add matrices in row and col structures & store in nrow structure
	for(i=0, nnz=0; i < nrows; i++){
		for(j=rowptr[i], k=colptr[i]; j < rowptr[i+1] && k < colptr[i+1];){
			if(rowind[j] == colind[k]){
				nrowval[nnz] = da_max(rowval[j], colval[k]);
				nrowind[nnz++] = rowind[j];
				j++;
				k++;
			} else if(rowind[j] < colind[k]){
				nrowval[nnz] = rowval[j];
				nrowind[nnz++] = rowind[j];
				j++;
			} else {
				nrowval[nnz] = colval[k];
				nrowind[nnz++] = colind[k];
				k++;
			}
		}
		for(; j < rowptr[i+1]; j++){
			nrowval[nnz] = rowval[j];
			nrowind[nnz++] = rowind[j];
		}
		for(; k < colptr[i+1]; k++){
			nrowval[nnz] = colval[k];
			nrowind[nnz++] = colind[k];
		}
		nrowptr[i+1] = nnz;
	}

	da_free((void**)&mat->rowptr, &mat->rowind, &mat->rowval,
			&mat->colptr, &mat->colind, &mat->colval, LTERM);
	mat->rowptr = nrowptr;
	mat->rowind = nrowind;
	mat->rowval = nrowval;

}




/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
int da_log2(idx_t a)
{
  ssize_t i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}


/*************************************************************************
* This function checks if the argument is a power of 2
**************************************************************************/
int da_ispow2(idx_t a)
{
  return (a == (1<<da_log2(a)));
}


/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
float da_flog2(float a)
{
  return log(a)/log(2.0);
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
double da_dlog2(double a)
{
  return log(a)/log(2.0);
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
val_t da_vlog2(val_t a)
{
  return log(a)/log(2.0);
}




#define COMPERRPRINT(ar, ac, av, br, bc, bv, rc, fr) \
	printf("%sa[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_VALTYPE \
	"] != b[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_VALTYPE "]", fr++?", ":"", rc?(ar):(ac), rc?(ac):(ar), av, rc?(br):(bc), rc?(bc):(br), bv);

#define COMPERRPRINT_A(ar, ac, bv, rc, fr) \
    printf("%s!a[" PRNT_IDXTYPE "," PRNT_IDXTYPE ",(" PRNT_VALTYPE ")]", fr++?", ":"", rc?(ar):(ac), rc?(ac):(ar), bv);

#define COMPERRPRINT_B(br, bc, av, rc, fr) \
	printf("%s!b[" PRNT_IDXTYPE "," PRNT_IDXTYPE ",(" PRNT_VALTYPE ")]", fr++?", ":"", rc?(br):(bc), rc?(bc):(br), av);

#define COMPERRPRINT_VAL(ar, ac, av, bv, rc, fr) \
    printf("%s![" PRNT_IDXTYPE "," PRNT_IDXTYPE ",(" PRNT_VALTYPE ", " PRNT_VALTYPE ")]", fr++?", ":"", rc?(ar):(ac), rc?(ac):(ar), av, bv);


/**
 * Compare two csr matrices and print out differences
 * 	\param doc1 first matrix to compare
 * 	\param doc2 second matrix to compare
 * 	\param eps Float max delta for value comparison
 * 	\param compVals Whether values should be compared
 */
void da_csrCompare(da_csr_t* doc1, da_csr_t* doc2, float eps, char compInds, char compVals){
	ssize_t j, k, ndiff = 0;
	idx_t i, l, fr;
	da_csr_t *a = NULL, *b = NULL;
	ptr_t *ptr1, *ptr2;
	idx_t *ind1, *ind2;
	val_t *val1, *val2;
	char rc;

	ASSERT((doc1->rowptr && doc2->rowptr) || (doc1->colptr && doc2->colptr));

	a = da_csr_Copy(doc1);
	b = da_csr_Copy(doc2);
	if(compInds){
        da_csr_SortIndices(a, DA_ROW);
        da_csr_SortIndices(b, DA_ROW);
	}

	if(a->rowptr){
		rc = 1;
		ptr1 = a->rowptr;
		ind1 = a->rowind;
		val1 = a->rowval;
		ptr2 = b->rowptr;
		ind2 = b->rowind;
		val2 = b->rowval;
	} else {
		rc = 0;
		ptr1 = a->rowptr;
		ind1 = a->rowind;
		val1 = a->rowval;
		ptr2 = b->rowptr;
		ind2 = b->rowind;
		val2 = b->rowval;
	}
	if(!val1) compVals = 0;

	if((a->nrows != b->nrows) || (a->ncols != b->ncols) || (ptr1[a->nrows] != ptr2[b->nrows]))
		printf("Matrix stats differ: A[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE
				"] != B[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_PTRTYPE "].\n",
				a->nrows, a->ncols, ptr1[a->nrows],
				b->nrows, b->ncols, ptr2[b->nrows]);
	printf("Differences: \n");
	if(compVals && !compInds){
        for(i=0; i < a->nrows && i < b ->nrows; i++){
            fr=0;
            for(j=ptr1[i], k=ptr2[i]; j < ptr1[i+1] && k < ptr2[i+1]; ++j, ++k){
                if(da_abs(val1[j] - val2[k]) > eps){
                    printf("%s[%d %zu %f %f]", fr++?", ":"", i, j-ptr1[i]+1, val1[j], val2[k]);
                    ndiff++;
                }
            }
            for( ; j < ptr1[i+1]; j++ ){
                printf("%s!b[%d %zu %f]", fr++?", ":"", i, j-ptr1[i]+1, val1[j]);
                ndiff++;
            }
            for( ; k < ptr2[i+1]; k++ ){
                printf("%s!a[%d %zu %f]", fr++?", ":"", i, k-ptr2[i]+1, val2[k]);
                ndiff++;
            }
            if(fr) printf("\n");
        }
        for(l=i; i < a->nrows; i++){
            for(fr=0, j=ptr1[i]; j < ptr1[i+1]; j++ ){
                printf("%s!b[%d %zu %f]", fr++?", ":"", i, j-ptr1[i]+1, val1[j]);
                ndiff++;
            }
            if(fr) printf("\n");
        }
        i=l;
        for( ; i < b->nrows; i++){
            for(fr=0, k=ptr2[i]; k < ptr2[i+1]; k++ ){
                printf("%s!a[%d %zu %f]", fr++?", ":"", i, k-ptr2[i]+1, val2[k]);
                ndiff++;
            }
            if(fr) printf("\n");
        }

	} else {

        for(i=0; i < a->nrows && i < b ->nrows; i++){
            fr=0;
            for(j=ptr1[i], k=ptr2[i]; j < ptr1[i+1] && k < ptr2[i+1]; ){
                if(ind1[j] == ind2[k]){
                    if(compVals && da_abs(val1[j] - val2[k]) > eps){
                        COMPERRPRINT(i+1, ind1[j]+1, val1[j], i+1, ind2[k]+1, val2[k], rc, fr);
                        ndiff++;
                    }
                    j++;
                    k++;
                } else if(ind1[j] > ind2[k]){
                    COMPERRPRINT_A(i+1, ind2[k]+1, val2[k], rc, fr);
                    k++;
                    ndiff++;
                } else {
                    COMPERRPRINT_B(i+1, ind1[j]+1, val1[j], rc, fr);
                    j++;
                    ndiff++;
                }
            }
            for( ; j < ptr1[i+1]; j++ ){
                COMPERRPRINT_B(i+1, ind1[j]+1, val1[j], rc, fr);
                ndiff++;
            }
            for( ; k < ptr2[i+1]; k++ ){
                COMPERRPRINT_A(i+1, ind2[k]+1, val2[k], rc, fr);
                ndiff++;
            }
            if(fr) printf("\n");
        }
        for(l=i; i < a->nrows; i++){
            for(fr=0, j=ptr1[i]; j < ptr1[i+1]; j++ ){
                COMPERRPRINT_B(i+1, ind1[j]+1, val1[j], rc, fr);
                ndiff++;
            }
            if(fr) printf("\n");
        }
        i=l;
        for( ; i < b->nrows; i++){
            for(fr=0, k=ptr2[i]; k < ptr2[i+1]; k++ ){
                COMPERRPRINT_A(i+1, ind2[k]+1, val2[k], rc, fr);
                ndiff++;
            }
            if(fr) printf("\n");
        }

	}

	printf("Overall, %zu differences were encountered between A and B.\n", ndiff);

	da_csr_FreeAll(&a, &b, LTERM);
}


/**
 * Compute 1/sqrt(x) on avg 4x faster than libm
 * http://eggroll.unbsj.ca/rsqrt/rsqrt.pdf
 */
float rsqrt32(float number) {
	uint32_t i;
	float x2, y;
	x2 = number * 0.5F;
	y = number;
	i = *(uint32_t *) &y;
	i = 0x5f375a86 - (i >> 1);
	y = *(float *) &i;
	y = y * (1.5F - (x2 * y * y));
	return y;
}

double rsqrt64(double number) {
	uint64_t i;
	double x2, y;
	x2 = number * 0.5;
	y = number;
	i = *(uint64_t *) &y;
	i = 0x5fe6eb50c7b537a9 - (i >> 1);
	y = *(double *) &i;
	y = y * (1.5 - (x2 * y * y));
	return y;
}


/**
 * Undo the matrix permutation
 */
void da_inversePermuteMatrix(da_csr_t **matP, idx_t* rowPerm, idx_t* colPerm)
{
	ssize_t i, j, k, orgI, orgJ, sz;
	ptr_t nnz;
	ptr_t *ptr, *nptr;
	idx_t *ind, *nind, *rperm, *cperm, *revperm = NULL, *itmp = NULL;
	val_t *val, *nval = NULL, *vtmp = NULL;
	da_csr_t *mat = *matP;

	if(mat->rowptr){

		sz      = mat->nrows;
		nnz     = mat->rowptr[sz];
		ptr     = mat->rowptr;
		ind     = mat->rowind;
		val     = mat->rowval;
		rperm   = rowPerm;
		cperm   = colPerm;
		nptr    = da_pmalloc(sz + 1, "da_inversePermuteMatrix: nptr");
		nind    = da_imalloc(nnz, "da_inversePermuteMatrix: nind");
		nval    = da_vmalloc(nnz, "da_inversePermuteMatrix: nval");
		revperm = da_imalloc(sz, "da_inversePermuteMatrix: revperm");
		nptr[0] = 0;

		for ( i=0; i<sz; i++ )
			revperm[rperm[i]] = i;

		for ( i=0; i<sz; i++ ){
			orgI = revperm[i];
			nptr[i+1] = nptr[i] + ( ptr[orgI+1] - ptr[orgI] );
			for ( j=nptr[i]; j<nptr[i+1]; j++ )
			{
				orgJ = ptr[orgI] + j - nptr[i];
				nind[j] = cperm[ind[orgJ]];
				if(val)
					nval[j] = val[orgJ];
			}
		}

		da_free((void**) &(*matP)->rowptr, &(*matP)->rowind, &(*matP)->rowval, LTERM);
		(*matP)->rowptr = nptr;
		(*matP)->rowind = nind;
		(*matP)->rowval = nval;

		da_csr_SortIndices(mat, DA_ROW);

		da_free ( (void**)&revperm, LTERM );
	}

	if(mat->colptr){

		sz      = mat->ncols;
		nnz     = mat->colptr[sz];
		ptr     = mat->colptr;
		ind     = mat->colind;
		val     = mat->colval;
		rperm   = colPerm;
		cperm   = rowPerm;
		nptr    = da_pmalloc(sz + 1, "da_inversePermuteMatrix: nptr");
		nind    = da_imalloc(nnz, "da_inversePermuteMatrix: nind");
		nval    = da_vmalloc(nnz, "da_inversePermuteMatrix: nval");
		revperm = da_imalloc(sz, "da_inversePermuteMatrix: revperm");
		nptr[0] = 0;

		for ( i=0; i<sz; i++ )
			revperm[rperm[i]] = i;

		for ( i=0; i<sz; i++ ){
			orgI = revperm[i];
			nptr[i+1] = nptr[i] + ( ptr[orgI+1] - ptr[orgI] );
			for ( j=nptr[i]; j<nptr[i+1]; j++ )
			{
				orgJ = ptr[orgI] + j - nptr[i];
				nind[j] = cperm[ind[orgJ]];
				if(val)
					nval[j] = val[orgJ];
			}
		}

		da_free((void**) &(*matP)->colptr, &(*matP)->colind, &(*matP)->colval, LTERM);
		(*matP)->colptr = nptr;
		(*matP)->colind = nind;
		(*matP)->colval = nval;

		da_csr_SortIndices((*matP), DA_COL);

		da_free ( (void**)&revperm, LTERM );
	}


}


/**
 * Undo the matrix permutation in place
 */
void da_inversePermuteRows(da_csr_t *mat, idx_t* rowPerm, idx_t* colPerm)
{
    ssize_t i, j, k, orgI, orgJ, sz;
    ptr_t nnz;
    ptr_t *ptr, *nptr;
    idx_t *ind, *nind, *rperm, *cperm, *revperm = NULL, *itmp = NULL;
    val_t *val, *nval = NULL, *vtmp = NULL;

    if(! rowPerm)
        return;

    if(!mat->rowptr)
        da_errexit("da_inversePermuteRows: row structure required!");

    sz      = mat->nrows;
    nnz     = mat->rowptr[sz];
    ptr     = mat->rowptr;
    ind     = mat->rowind;
    val     = mat->rowval;
    rperm   = rowPerm;
    cperm   = colPerm;
    nptr    = da_pmalloc(sz + 1, "da_inversePermuteMatrix: nptr");
    nind    = da_imalloc(nnz, "da_inversePermuteMatrix: nind");
    nval    = da_vmalloc(nnz, "da_inversePermuteMatrix: nval");
    revperm = da_imalloc(sz, "da_inversePermuteMatrix: revperm");
    nptr[0] = 0;

    for ( i=0; i<sz; i++ )
        revperm[rperm[i]] = i;

    for ( i=0; i<sz; i++ ){
        orgI = revperm[i];
        nptr[i+1] = nptr[i] + ( ptr[orgI+1] - ptr[orgI] );
        for ( j=nptr[i]; j<nptr[i+1]; j++ )
        {
            orgJ = ptr[orgI] + j - nptr[i];
            nind[j] = cperm[ind[orgJ]];
            if(val)
                nval[j] = val[orgJ];
        }
    }

    da_pcopy(sz+1, nptr, ptr);
    da_icopy(nnz, nind, ind);
    if(val)
        da_vcopy(nnz, nval, val);

    da_csr_SortIndices(mat, DA_ROW);

    da_free ( (void**)&revperm, &nptr, &nind, &nval, LTERM );

}


/**
 * Print compile choices for L2AP
 */
void printCompileChoices(){
	printf("types of pruning: ");
	int i=0;
#ifdef L2PS
	if(i++)printf("& "); printf("l2ps ");
#elif defined(LENPS)
	if(i++)printf("& "); printf("lenps ");
#endif
#if defined(L2PS) && defined(LENPS)
da_errexit("Only one of L2PS and LENPS can be used for index reduction.");
#endif

#ifdef SZ1
	if(i++)printf("& "); printf("sz1 ");
#elif defined(SZ3)
	if(i++)printf("& "); printf("sz3 ");
#endif

#ifdef RS3
	if(i++)printf("& "); printf("rs3 ");
#elif defined(RS1)
	if(i++)printf("& "); printf("rs1 ");
#endif
#ifdef RS2
	if(i++)printf("& "); printf("rs2 ");
#elif defined(RS4)
	if(i++)printf("& "); printf("rs4 ");
#endif
#if defined(RS2) && defined(RS4)
da_errexit("Only one of RS2 and RS4 can be used for candidate generation pruning.");
#endif

#ifdef L2CG
	if(i++)printf("& "); printf("l2cg ");
#elif defined(LENCG)
	if(i++)printf("& "); printf("lencg ");
#endif
#if defined(L2CG) && defined(LENCG)
da_errexit("Only one of L2CG and LENCG can be used for candidate generation pruning.");
#endif

#ifdef PSCV
	if(i++)printf("& "); printf("pscv ");
#endif

#ifdef DP1
	if(i++)printf("& "); printf("dp1 ");
#endif
#ifdef DP2
	if(i++)printf("& "); printf("dp2 ");
#endif
#ifdef DP3
	if(i++)printf("& "); printf("dp3 ");
#endif
#ifdef DP4
	if(i++)printf("& "); printf("dp4 ");
#endif
#ifdef DP5
	if(i++)printf("& "); printf("dp5 ");
#endif
#ifdef DP6
	if(i++)printf("& "); printf("dp6 ");
#endif
#ifdef DP7
	if(i++)printf("& "); printf("dp7 ");
#endif
#ifdef DP8
	if(i++)printf("& "); printf("dp8 ");
#endif

#ifdef L2CV
	if(i++)printf("& "); printf("l2cv ");
#elif defined(LENCV)
	if(i++)printf("& "); printf("lencv ");
#endif
#if defined(L2CV) && defined(LENCV)
da_errexit("Only one of L2CV and LENCV can be used for candidate verification pruning.");
#endif


	printf("\n");

}


/**
 * Save KNNG to file or print to screen. Acceptable output formats are
 * DA_FMT_CLUTO, DA_FMT_CSR, DA_FMT_IJV
 */
void da_save_knng(params_t *params, da_iapq_t** knng,
        const idx_t nrows, const char ftype)
{
    idx_t i, j;
    size_t n;

    if(!params->fpout)
        return;

    if(ftype != DA_FMT_CLUTO && ftype != DA_FMT_CSR && ftype != DA_FMT_IJV)
        da_errexit("Invalid output format. Acceptable output formats are "
                "DA_FMT_CLUTO, DA_FMT_CSR, and DA_FMT_IJV");

    for(n=0, i=0; i < nrows; ++i)
        n += knng[i]->nnodes;
    params->nSimPairs = n;
    if(ftype == DA_FMT_CLUTO)
        fprintf(params->fpout, "%d %d %zu\n", nrows, nrows, n);

    if(ftype == DA_FMT_CLUTO || ftype == DA_FMT_CSR) {
        for(i=0; i < nrows; ++i){
            da_iakvsortd(knng[i]->nnodes, knng[i]->heap);
            for(j=0; j < knng[i]->nnodes; ++j ){
                fprintf(params->fpout, PRNT_IDXTYPE " " PRNT_ACCUMTYPE " " ,
                        knng[i]->heap[j].key+1, knng[i]->heap[j].val);
            }
            fprintf(params->fpout, "\n");
        }
    } else {
        for(i=0; i < nrows; ++i){
            da_iakvsortd(knng[i]->nnodes, knng[i]->heap);
            for(j=0; j < knng[i]->nnodes; ++j ){
                fprintf(params->fpout, PRNT_IDXTYPE " " PRNT_IDXTYPE " " PRNT_ACCUMTYPE "\n" ,
                        i+1, knng[i]->heap[j].key+1, knng[i]->heap[j].val);
            }
        }
    }
    fflush(params->fpout);
}


/**
 * Transfer data from knng stored as heaps to the rows structure of a CSR matrix
 */
da_csr_t* knng2csr(da_iapq_t** knng, const idx_t nrows)
{
    size_t nnz, i, j;
    ptr_t *ptr;
    idx_t *ind;
    val_t *val;
    da_csr_t *mat = NULL;

    for(nnz=0, i=0; i < nrows; ++i)
        nnz += knng[i]->nnodes;

    mat = da_csr_Create();
    mat->nrows = mat->ncols = nrows;
    ptr = mat->rowptr = da_pmalloc(nrows+1, "knng2csr: mat->rowptr");
    ind = mat->rowind = da_imalloc(nnz, "knng2csr: mat->rowind");
    val = mat->rowval = da_vmalloc(nnz, "knng2csr: mat->rowval");

    for(nnz=0, ptr[0]=0, i=0; i < nrows; ++i){
        for(j=0; j < knng[i]->nnodes; ++j ){
            ind[nnz] = knng[i]->heap[j].key;
            val[nnz++] = knng[i]->heap[j].val;
        }
        ptr[i+1] = nnz;
    }

    return mat;
}

/**
 * Verify results given pre-computed stored results
 * \param ngbrs Neighbors found for each row in the input matrix
 * \param vfile File containing true neighbors, in CSR (no header) format,
 *              sorted in decreasing value order in each row
 */
void verify_knng_results(da_csr_t *ngbrs, char* vfile, char print_errors)
{
    float v, lv, lv2;
    double recall, crecall;
    size_t i, j, cid, ln, sz, c, cc, n, nrows, nrows2, ncols, nnz, lnlen, err;
    idx_t progressInd, pct;
    char *line=NULL, *head, *tail;
    FILE *fpin;

    idx_t *ind = ngbrs->rowind;
    val_t *val = ngbrs->rowval;
    ptr_t *ptr = ngbrs->rowptr;
    nrows = ngbrs->nrows;
    ncols = ngbrs->ncols;

    val_t *row = da_vsmalloc(nrows, -1.0, "row");

    da_getfilestats(vfile, &nrows2, &nnz, NULL, NULL);

    if(nrows != nrows2)
        da_errexit("Num rows in result %zd does not match that in the verification file %zd.\n", nrows, nrows2);
    if (nnz%2 == 1)
        da_errexit("Error: The number of numbers %zd in the input file is not even.\n", nnz);
    nnz = nnz/2;
    fpin = da_fopen(vfile, "r", "da_csr_Read: fpin");

    /* compare results, one row at a time */
    recall = crecall = 0.0;
    da_progress_init(pct, progressInd, nrows);
    printf("Checking recall... ");
    for (n=0, i=0; i<nrows; ++i) {
        if(ptr[i+1] == ptr[i]){
            if (i % progressInd == 0 )
                da_progress_advance(pct);
            da_getline(&line, &lnlen, fpin);
            continue;
        }
        for(lv=FLT_MAX, j=ptr[i]; j < ptr[i+1]; ++j){
            row[ind[j]] = val[j];
            if(val[j] < lv)
                lv = val[j];
        }

        do {
            if (da_getline(&line, &lnlen, fpin) == -1)
                da_errexit("Premature end of input file: file while reading row %d\n", i);
        } while (line[0] == '%');

        head = line;
        tail = NULL;
        ln = 0;
        sz = da_min(ptr[i+1]-ptr[i], (ssize_t)ncols);
        c = cc = err = 0;
        lv2 = FLT_MAX;
        while(1){
            cid = (int)strtol(head, &tail, 0);
            if (tail == head)
                break;
            head = tail;
#ifdef __MSC__
            v = (float)strtod(head, &tail);
#else
            v = strtof(head, &tail);
#endif
            if (tail == head)
                da_errexit("Value could not be found for column! Row:%zd, col:%zd\n", i, cid);
            head = tail;
            if(v < lv2)
                lv2 = v;
            if(row[cid-1] > -1){
                c++;
                if(da_abs(row[cid-1] - v) < 1e-4){
                    cc++;
                } else if(print_errors > 0){  /* show neighbors we found who's sim may be incorrectly computed */
                    printf("[%zu %zu %f %f] ", i+1, cid, v, row[cid-1]);
                    err++;
                }
                row[cid-1] = 1;
            } else if(da_abs(lv - v) < 1e-4){
                cc++;
                if(print_errors > 2){ /* show neighbors we did not find within the min values */
                    printf("[%zu *%zu %f] ", i+1, cid, v);
                    err++;
                }
            } else if(print_errors > 1){ /* show neighbors we did not find */
                printf("[%zu -%zu %f] ", i+1, cid, v);
                err++;
            }

            ln++;
            if(ln == sz)
                break;
        }
        if(ln > 0){
            recall += (double)c/(double)ln;
            crecall += (double)cc/(double)ln;
            n++;
        }
        for(j=ptr[i]; j < ptr[i+1]; ++j){
            if(print_errors > 2 && row[ind[j]] != 1 && da_abs(lv - row[ind[j]]) > 1e-4 ){  /* show extra neighbors we reported that are not in the true neighborhood */
                printf("[%zu +%d %f] ", i+1, ind[j]+1, row[ind[j]]);
                err++;
            }
            row[ind[j]] = -1;
        }
        if(print_errors > 0 && err){
            printf("min: %f %f %.5f\n", lv, lv2, da_abs(lv-lv2));
            fflush(stdout);
        }

        if (i % progressInd == 0 )
            da_progress_advance(pct);
    }
    da_progress_finalize(pct);
    da_fclose(fpin);
    da_free((void **)&line, &row, LTERM);

    printf("\nRecall: %.4f\n", recall/n);
    printf("Correct recall: %.4f\n", crecall/n);
}


void verify_knng_results2(da_iapq_t** knng, const idx_t nrows, char* vfile, char print_errors)
{
    da_csr_t *mat = knng2csr(knng, nrows);
    verify_knng_results(mat, vfile, print_errors);
    da_csr_Free(&mat);
}


/**
 * Verify results given pre-computed stored results
 * \param ngbrs1 True Neighbors found for each row in the input matrix
 * \param ngbrs2 Neighbors found for each row in the input matrix
 */
void verify_knng_results3(da_csr_t *ngbrs1, da_csr_t *ngbrs2, idx_t nsz, char print_errors)
{
    float v, lv, lv2;
    double recall, crecall;
    size_t i, j, k, cid, ln, sz, c, cc, n, nrows, nrows2, ncols, nnz, lnlen, err;
    idx_t progressInd, pct;
    char *line=NULL, *head, *tail;
    FILE *fpin;

    idx_t *ind = ngbrs1->rowind, *ind2 = ngbrs2->rowind;
    val_t *val = ngbrs1->rowval, *val2 = ngbrs2->rowval;
    ptr_t *ptr = ngbrs1->rowptr, *ptr2 = ngbrs2->rowptr;
    nrows = ngbrs1->nrows;

    if(nrows != ngbrs2->nrows)
        da_errexit("Num rows in result %xd does not match that in the verification file %zd.\n",
                nrows, ngbrs2->nrows);

    val_t *row = da_vsmalloc(nrows, -1.0, "row");

    /* compare results, one row at a time */
    recall = crecall = 0.0;
    da_progress_init(pct, progressInd, nrows);
    printf("Checking recall... ");
    for (n=0, i=0; i<nrows; ++i) {
        if(ptr[i+1] == ptr[i]){
            if (i % progressInd == 0 )
                da_progress_advance(pct);
            continue;
        }
        for(lv=FLT_MAX, j=ptr[i]; j < ptr[i+1]; ++j){
            row[ind[j]] = val[j];
            if(val[j] < lv)
                lv = val[j];
        }

        c = cc = err = 0;
        lv2 = FLT_MAX;
        ln = da_min(nsz, ptr2[i+1]-ptr2[i]);
        for(k=0, j=ptr2[i]; j < ptr2[i+1] && k < nsz; ++j, ++k){
            cid = ind2[j];
            v = val2[j];
            if(row[cid] > -1){
                c++;
                if(da_abs(row[cid] - v) < 1e-4){
                    cc++;
                } else if(print_errors > 0){  /* show neighbors we found who's sim may be incorrectly computed */
                    printf("[%zu %zu %f %f] ", i+1, cid+1, v, row[cid]);
                    err++;
                }
                row[cid] = 1;
            } else if(da_abs(lv - v) < 1e-4){
                cc++;
                if(print_errors > 1){ /* show neighbors we did not find within the min values */
                    printf("[%zu *%zu %f] ", i+1, cid+1, v);
                    err++;
                }
            } else if(print_errors > 1){ /* show neighbors we did not find */
                printf("[%zu -%zu %f] ", i+1, cid+1, v);
                err++;
            }
            if(v < lv2)
                lv2 = v;
        }
        if(ln > 0){
            recall += (double)c/(double)ln;
            crecall += (double)cc/(double)ln;
            n++;
        }
        for(k-=err, j=ptr[i]; j < ptr[i+1]; ++j){
            if(print_errors > 2 && k < nsz && row[ind[j]] != 1){  /* show extra neighbors we reported that are not in the true neighborhood */
                printf("[%zu +%d %f] ", i+1, ind[j]+1, row[ind[j]]);
                err++;
                k++;
            }
            row[ind[j]] = -1;
        }
        if(print_errors > 0 && err){
            printf("min: %f %f %.5f\n", lv, lv2, da_abs(lv-lv2));
            fflush(stdout);
        }

        if (i % progressInd == 0 )
            da_progress_advance(pct);
    }
    da_progress_finalize(pct);
    da_free((void **)&row, LTERM);

    printf("\nRecall: %.4f\n", recall/n);
    printf("Correct recall: %.4f\n", crecall/n);
}



void da_free_sims(da_sims_t** s){
    if((*s)->sims)
        da_free((void**)&(*s)->sims, LTERM);
    if((*s)->ptr)
        da_free((void**)&(*s)->ptr, LTERM);
    da_free((void**)s, LTERM);
}
