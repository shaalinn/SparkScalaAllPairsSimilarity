/*!
 \file  idxjoin.c
 \brief This file contains idxJoin related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
idx_t da_GetSimilarSmallerRows(da_csr_t *mat, idx_t rid, char noSelfSim,
	idx_t nqterms, idx_t *qind, val_t *qval, char simtype, idx_t nsim,
	float minsim, da_ivkv_t *hits, idx_t *i_marker, da_ivkv_t *i_cand, idx_t *ncands);

/**
 * Main entry point to IdxJoin.
 */
void ijFindNeighbors(params_t *params){

	ssize_t i, ii, c, j, k, nneighbs;
	idx_t nrows, ncands, progressInd;
	ptr_t *rowptr;
	idx_t *rowind, *marker = NULL;
	float simT;
	val_t *rowval;
	da_ivkv_t *hits = NULL, *cand = NULL;
	da_csr_t *docs, *neighbors = NULL;

	docs    = params->docs;
	nrows   = docs->nrows;  // num rows
	rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
	rowind  = docs->rowind; // colid for each nnz
	rowval  = docs->rowval; // val for each nnz
	simT    = params->simT;

	progressInd = nrows/20;
	ncands  = 0; // number of considered candidates

	// allocate memory
	gk_startwctimer(params->timer_5); // memory allocation time
    if(params->sim == DA_SIM_TAN && params->norm > 0){
        printf("Warning: setting -norm=0 since -sim=tan.\n");
        params->norm = 0;
    }
	simSearchSetup(params);
	if(params->sim == DA_SIM_TAN || params->sim == DA_SIM_JAC){
	    da_csr_ComputeSquaredNorms(docs, DA_ROW);
	}
	hits = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "apFindNeighbors: hits");
	cand = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "apFindNeighbors: cand");
	marker = da_ismalloc(nrows, -1, "apFindNeighbors: marker");
	gk_stopwctimer(params->timer_5); // memory allocation time

	gk_startwctimer(params->timer_7); // indexing time
	da_csr_CreateIndex(docs, DA_COL);
	gk_stopwctimer(params->timer_7); // indexing time

	if(params->nim)
		neighbors = params->neighbors;

	// for each x \in V do
	gk_startwctimer(params->timer_3); // find neighbors time
	if(params->verbosity > 0)
		printf("Progress Indicator (no. of records): ");
	for(i=0; i < nrows; i++){

		k = da_GetSimilarSmallerRows(docs, i, 1, rowptr[i+1]-rowptr[i],
				rowind + rowptr[i], rowval + rowptr[i], params->sim,
				nrows, simT, hits, marker, cand, &ncands);

		params->nCandidates += ncands;
		params->nDotProducts += ncands;

		//transfer candidates
		params->nSimPairs += k;
		if(params->nim){
			// set current num neighbors for this row's search
			neighbors->rowptr[i+1] = neighbors->rowptr[i];
			// grow neighborhood matrix if necessary
			if(neighbors->rowptr[i] + k > params->nghnnz){
				params->nghnnz += gk_max(k, params->nghnnz/2);
				da_csr_Grow(neighbors, params->nghnnz);
				params->nghinccnt++;
			}
		}

		for(j=0; j < k; j++)
			da_storeSim(params, i, hits[j].key, hits[j].val);

		if ( params->verbosity > 0 && i % progressInd == 0 ) {
			printf("%d%%..", (int) ceil(100*i/(float)nrows));
			fflush(stdout);
		}
	}
	if(params->verbosity > 0)
		printf("100%%\n");
	gk_stopwctimer(params->timer_3); // find neighbors time

	simSearchFinalize(params);
	gk_free((void**)&hits, &cand, &marker, LTERM);
}



/**
 * Find similar rows in the matrix -  this version of the function reports
 * the number of candidates/dot products that were considered in the search.
 */
idx_t da_GetSimilarSmallerRows(da_csr_t *mat, idx_t rid, char noSelfSim,
		idx_t nqterms, idx_t *qind, val_t *qval, char simtype, idx_t nsim,
        float minsim, da_ivkv_t *hits, idx_t *i_marker, da_ivkv_t *i_cand, idx_t *ncands)
{
	ssize_t i, ii, j, k;
	idx_t nrows, ncols, ncand;
	ptr_t *colptr;
	idx_t *colind, *marker;
	val_t *colval, *rnorms, *rsums, mysum;
	da_ivkv_t *cand;

	if (nqterms == 0)
		return 0;

	nrows  = mat->nrows;
	ncols  = mat->ncols;
    rnorms = mat->rnorms;
	colptr = mat->colptr;
	colind = mat->colind;
	colval = mat->colval;

	marker = (i_marker ? i_marker : da_ismalloc(nrows, -1, "da_csr_GetSimilarSmallerRows: marker"));
	cand   = (i_cand   ? i_cand   : da_ivkvmalloc(nrows, "da_csr_GetSimilarSmallerRows: cand"));

	switch (simtype) {
	case DA_SIM_COS:
		for (ncand=0, ii=0; ii<nqterms; ii++) {
			i = qind[ii];
			if (i < ncols) {
				for (j=colptr[i]; j<colptr[i+1]; j++) {
					k = colind[j];
					if(k > rid || (noSelfSim && k == rid))
						break;
					if (marker[k] == -1) {
						cand[ncand].key = k;
						cand[ncand].val = 0;
						marker[k]       = ncand++;
					}
					cand[marker[k]].val += colval[j]*qval[ii];
				}
			}
		}
		break;

    case DA_SIM_JAC:
    case DA_SIM_TAN:
		for (ncand=0, ii=0; ii<nqterms; ii++) {
			i = qind[ii];
			if (i < ncols) {
				for (j=colptr[i]; j<colptr[i+1]; j++) {
					k = colind[j];
					if(k > rid || (noSelfSim && k == rid))
						break;
					if (marker[k] == -1) {
						cand[ncand].key = k;
						cand[ncand].val = 0;
						marker[k]       = ncand++;
					}
					cand[marker[k]].val += colval[j]*qval[ii];
				}
			}
		}

		for (i=0; i<ncand; i++)
			cand[i].val /= (rnorms[cand[i].key]+rnorms[rid]-cand[i].val);
		break;


	default:
		gk_errexit(SIGERR, "Unknown similarity measure %d\n", simtype);
		return -1;
	}

	*ncands = ncand;

	/* go and prune the hits that are bellow minsim */
	for (j=0, i=0; i<ncand; i++) {
		marker[cand[i].key] = -1;
		if (cand[i].val >= minsim)
			cand[j++] = cand[i];
	}
	ncand = j;

	if (nsim == -1 || nsim >= ncand) {
		nsim = ncand;
	}
	else {
		nsim = gk_min(nsim, ncand);
		da_ivkvkselectd(ncand, nsim, cand);
		da_ivkvsortd(nsim, cand);
	}

	da_ivkvcopy(nsim, cand, hits);

	if (i_marker == NULL)
		gk_free((void **)&marker, LTERM);
	if (i_cand == NULL)
		gk_free((void **)&cand, LTERM);

	return nsim;
}


