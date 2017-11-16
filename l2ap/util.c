/*!
 \file  util.c
 \brief This file contains utilities for the application

 \author David C. Anastasiu
 */

#include "includes.h"

/*************************************************************************
* This function returns the key of a particular StringMap ID
**************************************************************************/
char* da_getStringKey(const gk_StringMap_t *strmap, const char id)
{
  int i;

  for (i=0; strmap[i].name; i++) {
    if (strmap[i].id == id)
      return strmap[i].name;
  }

  return NULL;
}


/*************************************************************************
* This function returns the ID of a particular string based on the
* supplied StringMap array
**************************************************************************/
int da_getStringID(const gk_StringMap_t *strmap, char *key)
{
  int i;

  for (i=0; strmap[i].name; i++) {
    if (gk_strcasecmp(key, strmap[i].name))
      return strmap[i].id;
  }

  return -1;
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
	} else if(gk_fexists(file)){ // assume some sort of CSR. Can we guess?
		gk_getfilestats(file, NULL, &nnz, NULL, NULL);
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
		fpout = gk_fopen(filename, "w", "da_vWriteVector: fpout");
	else
		fpout = stdout;

	for(i=0; i < size; i++)
		fprintf(fpout, PRNT_VALTYPE "%s", vec[i], separator);

	if (filename)
		gk_fclose(fpout);
}

void da_vWriteMatrix(char* filename, val_t** mat, idx_t nrows, idx_t ncols){
	idx_t i, j;
	FILE* fpout;

	if (filename)
		fpout = gk_fopen(filename, "w", "da_vWriteMatrix: fpout");
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
		gk_fclose(fpout);
}

void da_dWriteMatrix(char* filename, double** mat, idx_t nrows, idx_t ncols){
	idx_t i, j;
	FILE* fpout;

	if (filename)
		fpout = gk_fopen(filename, "w", "da_vWriteMatrix: fpout");
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
		gk_fclose(fpout);
}

void da_iWriteVector(char* filename, idx_t* vec, idx_t size, char *separator){
	idx_t i;
	FILE* fpout;

	if (filename)
		fpout = gk_fopen(filename, "w", "da_iWriteVector: fpout");
	else
		fpout = stdout;

	for(i=0; i < size; i++)
		fprintf(fpout, PRNT_IDXTYPE "%s", vec[i], separator);

	if (filename)
		gk_fclose(fpout);
}


void da_pWriteVector(char* filename, idx_t* vec, idx_t size, char *separator){
	idx_t i;
	FILE* fpout;

	if (filename)
		fpout = gk_fopen(filename, "w", "da_iWriteVector: fpout");
	else
		fpout = stdout;

	for(i=0; i < size; i++)
		fprintf(fpout, PRNT_IDXTYPE "%s", vec[i], separator);

	if (filename)
		gk_fclose(fpout);
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

	if(params->nim){
		nneighbs = params->neighbors->rowptr[i+1];
		params->neighbors->rowind[nneighbs] = j;
		params->neighbors->rowval[nneighbs] = v;
		params->neighbors->rowptr[i+1]++;
	} else {

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
			if(fpout)
                for(k=0; k < cache->nsims; k++){
                    fprintf(fpout, PRNT_IDXTYPE" "PRNT_IDXTYPE" %.8f\n"PRNT_IDXTYPE" "PRNT_IDXTYPE" %.8f\n", (int32_t) sims[k].i + zidx,
                            (int32_t) sims[k].j + zidx, (float) sims[k].val, (int32_t) sims[k].j + zidx,
                            (int32_t) sims[k].i + zidx, (float) sims[k].val);
                }
			cache->nsims = 0;
		}

	}
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
 * stopped with gk_startwctimer(t) and gk_stopwctimer(t) respectively.
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
 * gk_startwctimer(t) and gk_stopwctimer(t) respectively.

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
				nrowval[nnz] = gk_max(rowval[j], colval[k]);
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

	gk_free((void**)&mat->rowptr, &mat->rowind, &mat->rowval,
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

/**
 * Compare two csr matrices and print out differences
 * 	\param doc1 first matrix to compare
 * 	\param doc2 second matrix to compare
 * 	\param eps Float max delta for value comparison
 * 	\param compVals Whether values should be compared
 */
void da_csrCompare(da_csr_t* doc1, da_csr_t* doc2, float eps, char compVals){
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
	da_csr_SortIndices(a, DA_ROW);
	da_csr_SortIndices(b, DA_ROW);

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
	for(i=0; i < a->nrows && i < b ->nrows; i++){
		fr=0;
		for(j=ptr1[i], k=ptr2[i]; j < ptr1[i+1] && k < ptr2[i+1]; ){
			if(ind1[j] == ind2[k]){
				if(compVals && gk_abs(val1[j] - val2[k]) > eps){
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

		gk_free((void**) &(*matP)->rowptr, &(*matP)->rowind, &(*matP)->rowval, LTERM);
		(*matP)->rowptr = nptr;
		(*matP)->rowind = nind;
		(*matP)->rowval = nval;

		da_csr_SortIndices(mat, DA_ROW);

		gk_free ( (void**)&revperm, LTERM );
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

		gk_free((void**) &(*matP)->colptr, &(*matP)->colind, &(*matP)->colval, LTERM);
		(*matP)->colptr = nptr;
		(*matP)->colind = nind;
		(*matP)->colval = nval;

		da_csr_SortIndices((*matP), DA_COL);

		gk_free ( (void**)&revperm, LTERM );
	}


}


void printCompileChoices(){
	printf("types of pruning: ");
	int i=0;
#ifdef L2PS
	if(i++)printf("& "); printf("l2ps ");
#elif defined(LENPS)
	if(i++)printf("& "); printf("lenps ");
#endif
#if defined(L2PS) && defined(LENPS)
gk_errexit(SIGERR, "Only one of L2PS and LENPS can be used for index reduction.");
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
gk_errexit(SIGERR, "Only one of RS2 and RS4 can be used for candidate generation pruning.");
#endif

#ifdef L2CG
	if(i++)printf("& "); printf("l2cg ");
#elif defined(LENCG)
	if(i++)printf("& "); printf("lencg ");
#endif
#if defined(L2CG) && defined(LENCG)
gk_errexit(SIGERR, "Only one of L2CG and LENCG can be used for candidate generation pruning.");
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
gk_errexit(SIGERR, "Only one of L2CV and LENCV can be used for candidate verification pruning.");
#endif

#ifdef TL1
	if(i++)printf("& "); printf("tl1 ");
#endif
#ifdef TL2
	if(i++)printf("& "); printf("tl2 ");
#endif
#ifdef TRS
    if(i++)printf("& "); printf("trs ");
#endif



	printf("\n");

}
