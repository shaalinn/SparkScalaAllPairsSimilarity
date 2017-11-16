/*!
 \file  idxjoin.c
 \brief This file contains idxJoin related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations

/**
 * Main entry point to KnnL2ap.
 */
void kl2apFindNeighbors(params_t *params)
{

	ssize_t i, j, k, nneighbs, rs;
	idx_t nrows, rid, cid;
	val_t simT;
	ptr_t *nptr;
	idx_t *nind, *revperm, *rperm;
	val_t *nval;
	da_csr_t *docs, *neighbors=NULL;
	ptr_t *rloc=NULL;
    da_iapq_t **knng=NULL;
    char *querydocs;

	docs    = params->docs;
	nrows   = docs->nrows;  // num rows

    if(params->verbosity > 0)
        printf("Docs matrix: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, "
            PRNT_PTRTYPE " nnz\n", docs->nrows, docs->ncols, docs->rowptr[docs->nrows]);

	da_startwctimer(params->timer_3); // overall knn graph construction time

    // pre-process the data
    preProcessData(params);

    // compact the column space
    da_csr_CompactColumns(docs);

	// allocate memory
	da_startwctimer(params->timer_5); // memory allocation time
	params->rsizes = da_imalloc(docs->nrows, "l2apFindNeighbors: rsizes");
    params->csizes = da_ismalloc(docs->ncols, 0, "l2apFindNeighbors: csizes");
    params->rperm = da_imalloc(docs->nrows, "l2apFindNeighbors: rperm");
    params->cperm = da_imalloc(docs->ncols, "l2apFindNeighbors: cperm");
    params->rwgts = da_vsmalloc(docs->nrows, 0.0, "l2apFindNeighbors: rwgts");
    params->cwgts = da_vsmalloc(docs->ncols, 0.0, "l2apFindNeighbors: cwgts");
    revperm = da_imalloc(nrows, "da_inversePermuteMatrix: revperm");
	rloc   = da_psmalloc(nrows, -1, "l2knnFindNeighbors: rloc"); /* shared locator for all the priority queues */
	knng   = params->knng = (da_iapq_t**) da_malloc(nrows * sizeof(da_iapq_t*),
	            "knng"); /* KNNG result as set of priority queues */
    for (i=0; i < nrows; ++i)
        knng[i] = da_iapqCreateShared(params->k, nrows, rloc);
    params->knng = knng;
	da_stopwctimer(params->timer_5); // memory allocation time

    l2apReorderDocs(params, &docs, params->rsizes, params->csizes,
            params->rwgts, params->cwgts, params->rperm, params->cperm);
    rperm = params->rperm;
    for ( i=0; i<nrows; i++ )
        revperm[rperm[i]] = i;

	simT = 1.0;
	rs = nrows;

    if(params->cr){
        sprintf(params->crfname, "kl2ap-%s-%s-%d-%d-%d-%.3f.cr", da_getDataset(params),
            da_getStringKey(sim_options, params->sim), params->k,
            (int) params->scale, (int) params->norm, params->step);
        sprintf(params->crfname2, "%s0", params->crfname);
        if(da_fexists(params->crfname)){
            da_crRead(params, knng, &simT, NULL);
            // see which neighborhoods are done
            querydocs = params->querydocs = da_csmalloc(nrows, 1, "l2ap2FindNeighbors: querydocs");
            for(i=0; i < nrows; ++i){
                if( knng[i]->nnodes == params->k )
                    querydocs[revperm[i]] = 0;
            }

            rs = da_csum(nrows, querydocs, 1);
        }
    }

	while(simT > 0.0 && rs > 0){
	    simT -= params->step;
	    if(simT < 0)
	        simT = 0.0;
	    printf("simT %.3f\n", simT);

		params->indexSize = 0;
		params->simT = simT;

	    da_stopwctimer(params->timer_3); // overall knn graph construction time - l2apFindNeighbors starts and stops timers
		l2apFindNeighbors(params);
	    da_startwctimer(params->timer_3); // overall knn graph construction time

	    if(params->cr)
	        da_crWrite(params, knng, simT, 0);

		// see which neighborhoods are done
		querydocs = params->querydocs;
		for(i=0; i < nrows; ++i){
		    if( knng[i]->nnodes == params->k )
		        querydocs[revperm[i]] = 0;
		}

		rs = da_csum(nrows, querydocs, 1);
		printf("\nRemaining: %zu\n", rs);
	}

	da_stopwctimer(params->timer_3); // find neighbors time

	/* finalize search */

    /* verify results */
    if (params->vFile)
        verify_knng_results2(knng, nrows, params->vFile, params->verbosity);

    /* report number of similarity pairs found */
    for(i=0; i < nrows; ++i)
        params->nSimPairs += knng[i]->nnodes;

    /* save neighbor graph */
    da_save_knng(params, knng, nrows, DA_FMT_CLUTO);

    /* free memory */
    for (i=0; i < nrows; ++i) {
        knng[i]->locator = NULL;
        da_iapqDestroy(knng[i]);
    }
    da_free((void**) &knng, &rloc, &revperm, LTERM);
}
