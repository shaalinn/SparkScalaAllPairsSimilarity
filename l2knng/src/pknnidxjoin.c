/*!
 \file  idxjoin.c
 \brief This file contains parallel kIdxJoin related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
void ompFindNeighbors(params_t *params, da_csr_t *mat, da_csr_t *pmat,
        int qID, int nqrows,
         int dID, int ndrows,
         idx_t *nallhits, idx_t *nallcands, da_ivkv_t ** allhits,
         da_ivkv_t **atcand, da_ivkv_t **athits, idx_t **atmarker,
         char simtype);
#define NITERS          20

/**
 * Main entry point to ParallelKnnIdxJoin.
 */
void pkijFindNeighbors(params_t *params)
{

	ssize_t i, j, k, nneighbs;
	idx_t nrows, ncands, progressInd, pct, nqrows, ndrows, qID, dID;
	idx_t **marker=NULL, *nallhits=NULL, *nallcands=NULL, *csizes=NULL, *cperm=NULL;
	da_ivkv_t **hits=NULL, **cand=NULL, **allhits=NULL;
	da_csr_t *docs, *pmat;

	if(params->fpout && params->fmtWrite != DA_FMT_IJV)
	    params->nim = 1;

	docs    = params->docs;
	nrows   = docs->nrows;  // num rows
    params->nqrows = da_min(nrows, params->nqrows);
    params->ndrows = da_min(nrows, params->ndrows);

	if(params->fpout && params->fmtWrite != DA_FMT_CSR && params->fmtWrite != DA_FMT_IJV)
	    da_errexit("Mode pkij can only output data in CSR or IJV formats.");

	da_startwctimer(params->timer_3); // overall knn graph construction time
	/* pre-process the data */
	preProcessData(params);
    da_csr_CompactColumns(docs);

    if(params->permcol != PERM_NONE){
        csizes = da_ismalloc(docs->ncols, 0, "csizes");
        cperm  = da_imalloc(docs->ncols, "cperm");
    }
    if(params->permcol == PERM_L2AP){
        printf("colperm l2ap\n");
        l2apReorderCols(docs, csizes, cperm);
        da_free((void**)&cperm, &csizes, LTERM);
    } else if(params->permcol == PERM_L2KNN){
        printf("colperm l2knn\n");
        l2knnReorderCols(docs, csizes, cperm);
        da_free((void**)&cperm, &csizes, LTERM);
    }
    if(params->sim == DA_SIM_JAC)
        da_csr_ComputeSquaredNorms(docs, DA_ROW);

	/* allocate memory */
	da_startwctimer(params->timer_5); // memory allocation time
    da_AllocMatrix((void ***)&hits, sizeof(da_ivkv_t), params->nthreads, params->ndrows);
    da_AllocMatrix((void ***)&cand, sizeof(da_ivkv_t), params->nthreads, params->ndrows);
    marker = da_iAllocMatrix(params->nthreads, params->ndrows, -1, "pkijFindNeighbors: marker");
    da_AllocMatrix((void ***)&allhits, sizeof(da_ivkv_t), params->nqrows, 2*params->k);
    nallhits = da_imalloc(params->nqrows, "pkijFindNeighbors: nallhits");
    nallcands = da_imalloc(params->nqrows, "pkijFindNeighbors: nallcands");
	da_stopwctimer(params->timer_5); // memory allocation time
    da_stopwctimer(params->timer_3); // find neighbors time

	omp_set_num_threads(params->nthreads);

	// for each x \in V do
	if(params->verbosity > 0)
		printf("Progress Indicator: ");
	fflush(stdout);
    da_progress_init(pct, progressInd, nrows/params->nqrows);
	for (qID=0; qID < nrows; qID+=params->nqrows) {
	    da_startwctimer(params->timer_3); // overall knn graph construction time
	    nqrows = da_min(params->nqrows, nrows-qID);
	    da_iset(nqrows, 0, nallhits);

	    if (params->verbosity > 3)
	        printf("Working on query chunk: %7d, %4d\n", qID, nqrows);

	    /* find the neighbors of the chunk */
	    for (dID=0; dID < nrows; dID+=params->ndrows) {
	        ndrows = da_min(params->ndrows, nrows-dID);

	        /* create the sub-matrices */
	        pmat = da_csr_ExtractSubmatrix(docs, dID, ndrows);
	        ASSERT(pmat != NULL);
	        da_csr_CreateIndex(pmat, DA_COL);

	        if (params->verbosity > 4)
	            printf("  Working on db chunk: %7d, %4d, %4.2fMB\n", dID, ndrows,
	                    8.0*pmat->rowptr[pmat->nrows]/(1024*1024));

	        /* spawn the work threads */
	        ompFindNeighbors(params, docs, pmat, qID, nqrows, dID, ndrows,
	                nallhits, nallcands, allhits,
	                cand, hits, marker, params->sim);

	        params->nCandidates += da_isum(params->nqrows, nallcands, 1);

	        da_csr_Free(&pmat);
	    }
	    da_stopwctimer(params->timer_3); // find neighbors time

	    /* write the results in the file */
	    if (params->fpout) {
	        if (params->verbosity > 4)
	            printf(" writing... \n");
	        if(params->fmtWrite == DA_FMT_CSR)
                for (i=0; i<nqrows; i++) {
                    for (j=0; j<nallhits[i]; j++) {
                        fprintf(params->fpout, " %d %.7g", allhits[i][j].key+1, allhits[i][j].val);
                    }
                    fprintf(params->fpout, "\n");
                }
	        else
	            for (i=0; i<nqrows; i++) {
                    for (j=0; j<nallhits[i]; j++) {
                        fprintf(params->fpout, "%zu %d %.7g\n", qID+i+1, allhits[i][j].key+1, allhits[i][j].val);
                    }
                }
	        fflush(params->fpout);
	    }
        params->nSimPairs += da_isum(params->nqrows, nallhits, 1);

        if ( params->verbosity > 0 && (qID/params->nqrows) % progressInd == 0 )
            da_progress_advance(pct);

	}
	if(params->verbosity > 0){
        da_progress_finalize(pct);
	    printf("\n");
	}
	params->nDotProducts = params->nCandidates;

	/* free memory */
    da_FreeMatrix((void ***)&hits, params->nthreads, params->ndrows);
    da_FreeMatrix((void ***)&cand, params->nthreads, params->ndrows);
    da_FreeMatrix((void ***)&allhits, params->nqrows, 2*params->k);
    da_iFreeMatrix(&marker, params->nthreads, params->ndrows);
	da_free((void**)&nallhits, &nallcands, LTERM);
}




/*************************************************************************/
/*! Computes the neighbors of a set of rows against the documents in
    vault->pmat using OpenMP */
/**************************************************************************/
void ompFindNeighbors(params_t *params, da_csr_t *mat, da_csr_t *pmat,
        int qID, int nqrows,
         int dID, int ndrows,
         idx_t *nallhits, idx_t *nallcands, da_ivkv_t ** allhits,
         da_ivkv_t **atcand, da_ivkv_t **athits, idx_t **atmarker,
         char simtype)
{

    #pragma omp parallel
    {
        int i, j, k, l, ci, nchits, nhits, ncands, nnbrs, noldhits, tid;
        idx_t *marker;
        da_ivkv_t *cand, *hits, *newhits, *oldhits;

        tid = omp_get_thread_num();

        marker  = atmarker[tid];
        cand    = atcand[tid];
        hits   = athits[tid];

        #pragma omp for schedule(dynamic, NITERS)
        for (i=0; i < nqrows; i++) {
            /* compute the similarity */
            nhits = da_GetKSimilarRows(pmat, qID-dID+i, 1,
                    mat->rowptr[qID+i+1]-mat->rowptr[qID+i],
                    mat->rowind+mat->rowptr[qID+i],
                    mat->rowval+mat->rowptr[qID+i],
                    simtype, params->k, hits, marker, cand, &nallcands[i]);

            for (k=0; k<nhits; k++)
                hits[k].key += dID;

            /* merge with the current best neighbors */
            da_ivkvsortd(nhits, hits);
            nnbrs = da_min(params->k, nhits);

            newhits  = allhits[i];
            oldhits  = newhits + params->k;
            noldhits = nallhits[i];
            memcpy(oldhits, newhits, sizeof(da_ivkv_t)*noldhits);

            /* the two lists to be merged are (nnbrs, hits) and (noldhits, oldhits) */
            for (l=0, j=0, k=0; j < nnbrs && k < noldhits && l < params->k; l++) {
                if (hits[j].val >= oldhits[k].val) {
                    newhits[l] = hits[j];
                    j++;
                }
                else {
                    newhits[l] = oldhits[k];
                    k++;
                }
            }
            for (; j<nnbrs && l < params->k; l++, j++)
                newhits[l] = hits[j];
            for (; k<noldhits && l < params->k; l++, k++)
                newhits[l] = oldhits[k];
            nallhits[i] = l;

        }

    }

}


