/*!
 \file  knnidxjoin.c
 \brief This file contains kIdxJoin related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations

/**
 * Main entry point to KnnIdxJoin.
 */
void kijFindNeighbors(params_t *params)
{

	ssize_t i, j, k, nneighbs;
	size_t rid;
	idx_t nrows, ncands, progressInd, pct;
	ptr_t *rowptr;
	idx_t *rowind, *marker=NULL;
	val_t *rowval;
	da_ivkv_t *hits=NULL, *cand=NULL;
	idx_t *csizes=NULL, *cperm=NULL;
	da_csr_t *docs, *neighbors=NULL;

	if(params->fpout && params->fmtWrite != DA_FMT_IJV)
	    params->nim = 1;

	docs    = params->docs;
	nrows   = docs->nrows;  // num rows
	rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
	rowind  = docs->rowind; // colid for each nnz
	rowval  = docs->rowval; // val for each nnz

	params->nim = 1; /** temporary while running experiments **/

        da_progress_init_steps(pct, progressInd, nrows, 100);
	ncands  = 0; // number of considered candidates

	da_startwctimer(params->timer_3); // overall knn graph construction time

	// allocate memory
	da_startwctimer(params->timer_5); // memory allocation time
	simSearchSetup(params);
	hits = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "kijFindNeighbors: hits");
	cand = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "kijFindNeighbors: cand");
	marker = da_ismalloc(nrows, -1, "kijFindNeighbors: marker");
	if(params->permcol != PERM_NONE){
	    csizes = da_ismalloc(docs->ncols, 0, "csizes");
	    cperm  = da_imalloc(docs->ncols, "cperm");
	}
	da_stopwctimer(params->timer_5); // memory allocation time

	if(params->permcol == PERM_L2AP){
	    l2apReorderCols(docs, csizes, cperm);
	    da_free((void**)&cperm, &csizes, LTERM);
	} else if(params->permcol == PERM_L2KNN){
	    l2knnReorderCols(docs, csizes, cperm);
        da_free((void**)&cperm, &csizes, LTERM);
	}

	da_startwctimer(params->timer_7); // indexing time
	da_csr_CreateIndex(docs, DA_COL);
	da_stopwctimer(params->timer_7); // indexing time

	if(params->nim){
		neighbors = params->neighbors;
		params->nghnnz = nrows * params->k;
		da_csr_Grow(neighbors, params->nghnnz);
	}

	if(params->sim == DA_SIM_JAC)
	    da_csr_ComputeSquaredNorms(docs, DA_ROW);

	// for each x \in V do
	if(params->verbosity > 0)
		printf("Progress Indicator: ");
	i=0;

	if(params->cr){
	    sprintf(params->crfname, "kij-%s-%s-%d-%d-%d.cr", da_getDataset(params),
            da_getStringKey(sim_options, params->sim), params->k, (int) params->scale, (int) params->norm);
	    sprintf(params->crfname2, "%s0", params->crfname);
	    if(da_fexists(params->crfname)){
	        da_crRead(params, NULL, NULL, &rid);
	        i = rid;
            pct = 100*(i+1.0)/(float)nrows;
	    }
	}

	for(; i < nrows; i++){
		k = da_GetKSimilarRows(docs, i, 1, rowptr[i+1]-rowptr[i],
				rowind + rowptr[i], rowval + rowptr[i], params->sim,
				params->k, hits, marker, cand, &ncands);

		params->nCandidates += ncands;
		params->nDotProducts += ncands;

		// transfer candidates
		params->nSimPairs += k;
		if(params->nim)
			// set current num neighbors for this row's search
			neighbors->rowptr[i+1] = neighbors->rowptr[i];

		for(j=0; j < k; j++)
			da_storeSim(params, i, hits[j].key, hits[j].val);

		if ( params->verbosity > 0 && i % progressInd == 0 ){
		    if(params->cr)
                da_crWrite(params, NULL, 0, i);
            da_progress_advance_steps(pct, 100);
		}
	}
	if(params->verbosity > 0){
            da_progress_finalize_steps(pct, 100);
	    printf("\n");
	}
	if(!params->nim){
        // write out neighbors buffer
        da_storeSim(params, -1, -1, 0);
    }
	da_stopwctimer(params->timer_3); // find neighbors time

	/* finalize search */
    if(params->nim && params->fpout){
        if(params->verbosity > 0){
            printf("\nSorting neighbors in decreasing similarity value order...\n");
            fflush(stdout);
        }
        for(i=0; i < nrows; ++i){
            for(k=0, j=neighbors->rowptr[i]; j < neighbors->rowptr[i+1]; ++j, ++k){
                cand[k].key = neighbors->rowind[j];
                cand[k].val = neighbors->rowval[j];
            }
            da_ivkvsortd( neighbors->rowptr[i+1]-neighbors->rowptr[i], cand);
            for(k=0, j=neighbors->rowptr[i]; j < neighbors->rowptr[i+1]; ++j, ++k){
                neighbors->rowind[j] = cand[k].key;
                neighbors->rowval[j] = cand[k].val;
            }
        }

        if(params->verbosity > 0)
            printf("Writing neighborhood matrix to %s.\n", params->oFile);
        da_csr_Write(params->neighbors, params->oFile, params->fmtWrite, 1, 1);
    }
	da_free((void**)&hits, &cand, &marker, LTERM);
}

/**
 * Find knng of a subset of rows in the matrix
 * \param docs Matrix to find neighbors in
 * \param which List of row ids to find neighbors for
 * \param wsize Size of array which
 * \param k Number of neighbors to return
 * \param sim Which similarity to compute - DA_SIM_COS, DA_SIM_JAC, etc.
 * \param verbosity Whether to output progress and timing information to screen
 */
da_sims_t* cknngIdxJoin(da_csr_t *docs, int32_t* which, ssize_t wsize, int32_t k,
        char sim, char verbosity)
{
    ssize_t i, j, r, l, nnbrs;
    idx_t nrows, ncands, progressInd, pct;
    ptr_t *rowptr;
    idx_t *rowind, *marker = NULL;
    val_t *rowval;
    da_ivkv_t *hits = NULL, *cand = NULL;
    double timer = 0.0;
    da_sims_t *neighbors;


    nrows   = docs->nrows;  // num rows
    rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
    rowind  = docs->rowind; // colid for each nnz
    rowval  = docs->rowval; // val for each nnz

    da_progress_init(pct, progressInd, wsize);
    ncands  = 0; // number of considered candidates

    da_startwctimer(timer); // overall knn graph construction time

    // allocate memory
    hits = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "cknng: hits");
    cand = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "cknng: cand");
    marker = da_ismalloc(nrows, -1, "cknng: marker");
    neighbors = (da_sims_t*) da_malloc(sizeof(da_sims_t), "cknng: neighbors");
    neighbors->nsims = wsize;
    neighbors->size  = wsize * k;
    neighbors->ptr   = da_pmalloc(wsize+1, "cknng: neighbors->ptr");
    neighbors->sims  = da_smalloc(neighbors->size, "cknng: neighbors->ptr");
    neighbors->ptr[0] = 0;

    if(!docs->colptr)
        da_csr_CreateIndex(docs, DA_COL);

    if(sim == DA_SIM_JAC)
        da_csr_ComputeSquaredNorms(docs, DA_ROW);

    // for each x \in V do
    for(i=0; i < wsize; i++){
        r = which[i];
        nnbrs = da_GetKSimilarRows(docs, r, 1, rowptr[r+1]-rowptr[r],
                rowind + rowptr[r], rowval + rowptr[r], sim,
                k, hits, marker, cand, &ncands);

        // transfer candidates
        for(l=neighbors->ptr[i], j=0; j < nnbrs; ++j, ++l){
            neighbors->sims[l].i = r;
            neighbors->sims[l].j = hits[j].key;
            neighbors->sims[l].val = hits[j].val;
        }
        neighbors->ptr[i+1] = neighbors->ptr[i] + nnbrs;

        if ( verbosity > 0 && i % progressInd == 0 )
            da_progress_advance(pct);
    }
    da_stopwctimer(timer); // find neighbors time

    if(verbosity > 0){
        da_progress_finalize(pct);
        printf(" sz %zu - ", wsize);
        da_printTimerLong("",
            da_getwctimer(timer));
        fflush(stdout);
    }

    da_free((void**)&hits, &cand, &marker, LTERM);

    return neighbors;
}



/**
 * Find knng of a subset of rows in the matrix and store results directly in appropriate heaps in knng
 * \param docs Matrix to find neighbors in
 * \param which List of row ids to find neighbors for
 * \param wsize Size of array which
 * \param k Number of neighbors to return
 * \param sim Which similarity to compute - DA_SIM_COS, DA_SIM_JAC, etc.
 * \param verbosity Whether to output progress and timing information to screen
 */
void cknngIdxJoin2(da_iapq_t **knng, da_csr_t *docs, int32_t* which, ssize_t wsize, int32_t k,
        char sim, char verbosity)
{
    ssize_t i, j, r, l, nnbrs;
    idx_t nrows, ncands, progressInd, pct;
    ptr_t *rowptr;
    idx_t *rowind, *marker = NULL;
    val_t *rowval;
    da_ivkv_t *hits = NULL, *cand = NULL;
    double timer = 0.0;
    da_iapq_t *ngbr;

    nrows   = docs->nrows;  // num rows
    rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
    rowind  = docs->rowind; // colid for each nnz
    rowval  = docs->rowval; // val for each nnz


    da_progress_init(pct, progressInd, wsize);
    ncands  = 0; // number of considered candidates

    da_startwctimer(timer); // overall knn graph construction time

    // allocate memory
    hits = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "cknng: hits");
    cand = da_ivkvsmalloc(nrows, (da_ivkv_t) {0, 0.0}, "cknng: cand");
    marker = da_ismalloc(nrows, -1, "cknng: marker");

    if(!docs->colptr)
        da_csr_CreateIndex(docs, DA_COL);

    if(sim == DA_SIM_JAC)
        da_csr_ComputeSquaredNorms(docs, DA_ROW);

    // for each x \in V do
    for(i=0; i < wsize; i++){
        r = which[i];
        nnbrs = da_GetKSimilarRows(docs, r, 1, rowptr[r+1]-rowptr[r],
                rowind + rowptr[r], rowval + rowptr[r], sim,
                k, hits, marker, cand, &ncands);

        // transfer candidates
        ngbr = knng[r];
        ngbr->nnodes = nnbrs;
        for(l=0, j=0; j < nnbrs; ++j, ++l){
            ngbr->heap[l].key = hits[j].key;
            ngbr->heap[l].val = hits[j].val;
        }

        if ( verbosity > 0 && i % progressInd == 0 )
            da_progress_advance(pct);
    }
    da_stopwctimer(timer); // find neighbors time

    if(verbosity > 0){
        da_progress_finalize(pct);
        printf(" sz %zu - ", wsize);
        da_printTimerLong("",
            da_getwctimer(timer));
        fflush(stdout);
    }

    da_free((void**)&hits, &cand, &marker, LTERM);

}




/**
 * Find similar rows in the matrix -  this version of the function reports
 * the number of candidates/dot products that were considered in the search.
 * \param mat The CSR matrix we're searching in
 * \param rid Row we're looking for neighbors for
 * \param noSelfSim Do not return self-similarity as one of the results
 * \param nqterms length of query sparse vector
 * \param qind Column indices of the query sparse vector
 * \param qval Values in the query sparse vector
 * \param simtype Which similarity to compute (DA_SIM_COS or DA_SIM_JAC)
 * \param nsim Number of similar pairs to get (-1 to get all)
 * \param hits Array or length mat->nrows to hold values for possible matches and result
 * \param i_marker Optional marker array of length mat->nrows to mark candidates
 * \param i_cand Optional key-value array of length mat->nrows to store and sort candidates
 * \param ncands Reference to int variable to hold number of candidates
 *
 * \return Number of similar pairs found
 */
idx_t da_GetKSimilarRows(da_csr_t *mat, idx_t rid, char noSelfSim,
		idx_t nqterms, idx_t *qind, val_t *qval, char simtype, idx_t nsim,
        da_ivkv_t *hits, idx_t *i_marker, da_ivkv_t *i_cand, idx_t *ncands)
{
	ssize_t i, ii, j, k;
	idx_t nrows, ncols, ncand;
	ptr_t *colptr;
	idx_t *colind, *marker;
	val_t *colval, *rnorms, mynorm, *rsums, mysum;
	da_ivkv_t *cand;

	if (nqterms == 0)
		return 0;

	nrows  = mat->nrows;
	ncols  = mat->ncols;
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
					if(noSelfSim && k == rid)
						continue;
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
		for (ncand=0, ii=0; ii<nqterms; ii++) {
			i = qind[ii];
			if (i < ncols) {
				for (j=colptr[i]; j<colptr[i+1]; j++) {
					k = colind[j];
					if(noSelfSim && k == rid)
						continue;
					if (marker[k] == -1) {
						cand[ncand].key = k;
						cand[ncand].val = 0;
						marker[k]       = ncand++;
					}
					cand[marker[k]].val += colval[j]*qval[ii];
				}
			}
		}

		rnorms = mat->rnorms;
		mynorm = da_vdot(nqterms, qval, 1, qval, 1);

		for (i=0; i<ncand; i++)
			cand[i].val = cand[i].val/(rnorms[cand[i].key]+mynorm-cand[i].val);
		break;


	default:
		da_errexit("Unknown similarity measure %d\n", simtype);
		return 0;
	}

	*ncands = ncand;

	/* clear markers */
	for (j=0, i=0; i<ncand; i++)
        marker[cand[i].key] = -1;

	if (nsim == -1 || nsim >= ncand) {
		nsim = ncand;
	}
	else {
		nsim = da_min(nsim, ncand);
		da_ivkvkselectd(ncand, nsim, cand);
	}
	da_ivkvcopy(nsim, cand, hits);

	if (i_marker == NULL)
		da_free((void **)&marker, LTERM);
	if (i_cand == NULL)
		da_free((void **)&cand, LTERM);

	return nsim;
}


