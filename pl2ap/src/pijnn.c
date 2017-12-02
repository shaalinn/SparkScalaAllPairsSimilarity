/*!
 \file  pijnn.c
 \brief This file contains parallel IdxJoin Nearest Neighbor Search functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
void da_ompijnnFindNeighbors (params_t * params, da_csr_t * mat,
        da_csr_t * pmat, int qID, int nqrows, int dID,
        int ndrows, idx_t * nallhits, idx_t * nallcands,
        da_ivkv_t ** allhits, da_ivkv_t ** atcand,
        da_ivkv_t ** athits, idx_t ** atmarker,
        char simtype);
void da_ijnnPreProcessData (params_t * params);




/**
 * Main entry point to Parallel IdxJoin Nearest Neighbor Search
 */
void pijnnFindNeighbors (params_t * params)
{

    ssize_t i, j, k, nneighbs;
    idx_t nrows, nrowsdb, ncands, progressInd, pct, nqrows, ndrows, qID, dID;
    idx_t **marker = NULL, *nallhits = NULL, *nallcands = NULL, *csizes =
            NULL, *cperm = NULL;
    da_ivkv_t **hits = NULL, **cand = NULL, **allhits = NULL;
    da_csr_t *docs = NULL, *docsdb = NULL, *pmat = NULL;

    if (params->fpout && params->fmtWrite != DA_FMT_IJV)
        params->nim = 1;

    if (params->fpout && params->fmtWrite != DA_FMT_CSR && params->fmtWrite != DA_FMT_IJV)
        da_errexit ("Mode pijnn can only output data in CSR or IJV formats.");
    printf ("fmtWrite: %s\n", da_getStringKey (fmt_options, params->fmtWrite));

    da_startwctimer (params->timer_3);	// overall sim graph construction time
    /* pre-process the data */
    da_ijnnPreProcessData (params);

    docs    = params->docs;
    docsdb  = params->docsdb;
    nrows   = docs->nrows;       // num query rows
    nrowsdb = docsdb->nrows;       // num db rows
    params->nqrows = da_min (nrows, params->nqrows);
    params->ndrows = da_min (nrowsdb, params->ndrows);

    if (params->sim == DA_SIM_JAC){
        da_csr_ComputeSquaredNorms (docs, DA_ROW);
        da_csr_ComputeSquaredNorms (docsdb, DA_ROW);
    }

    /* allocate memory */
    da_startwctimer (params->timer_5);	// memory allocation time
    da_AllocMatrix ((void ***) &hits, sizeof (da_ivkv_t), params->nthreads,
            params->ndrows);
    da_AllocMatrix ((void ***) &cand, sizeof (da_ivkv_t), params->nthreads,
            params->ndrows);
    marker = da_iAllocMatrix (params->nthreads, params->ndrows, -1,
                    "pijFindNeighbors: marker");
    da_AllocMatrix ((void ***) &allhits, sizeof (da_ivkv_t), params->nqrows,
            nrowsdb + params->ndrows);
    nallhits = da_imalloc (params->nqrows, "pijFindNeighbors: nallhits");
    nallcands = da_imalloc (params->nqrows, "pijFindNeighbors: nallcands");
    da_stopwctimer (params->timer_5);	// memory allocation time
    da_stopwctimer (params->timer_3);	// find neighbors time

    omp_set_num_threads (params->nthreads);

    // for each x \in V do
    if (params->verbosity > 0)
        printf ("Progress Indicator: ");
    fflush (stdout);
    da_progress_init (pct, progressInd, (nrows / params->nqrows));
    for (qID = 0; qID < nrows; qID += params->nqrows) {
        da_startwctimer (params->timer_3);	// similarity search time -- includes indexing, cg and cv times
        nqrows = da_min (params->nqrows, nrows - qID);
        da_iset (nqrows, 0, nallhits);

        if (params->verbosity > 3)
            printf ("Working on query chunk: %7d, %4d\n", qID, nqrows);

        /* find the neighbors in the chunk */
        for (dID = 0; dID < nrowsdb; dID += params->ndrows) {
            ndrows = da_min (params->ndrows, nrowsdb - dID);

            /* create the database sub-matrix */
            pmat = da_csr_ExtractSubmatrix (docsdb, dID, ndrows);
            ASSERT (pmat != NULL);
            da_csr_CreateIndex (pmat, DA_COL);

            if (params->verbosity > 4)
                printf ("  Working on db chunk: %7d, %4d, %4.2fMB\n", dID, ndrows,
                        8.0 * pmat->rowptr[pmat->nrows] / (1024 * 1024));

            /* spawn the work threads */
            da_ompijnnFindNeighbors (params, docs, pmat, qID, nqrows, dID, ndrows,
                    nallhits, nallcands, allhits,
                    cand, hits, marker, params->sim);

            params->nCandidates += da_isum (params->nqrows, nallcands, 1);

            da_csr_Free (&pmat);
        }
        da_stopwctimer (params->timer_3);	// find neighbors time

        /* write the results in the file */
        if (params->fpout) {
            if (params->verbosity > 4)
                printf (" writing... \n");
            if (params->fmtWrite == DA_FMT_CSR)
                for (i = 0; i < nqrows; i++) {
                    for (j = 0; j < nallhits[i]; j++) {
                        fprintf (params->fpout, " %d %.7g",
                                allhits[i][j].key + 1, allhits[i][j].val);
                    }
                    fprintf (params->fpout, "\n");
                } else
                    for (i = 0; i < nqrows; i++) {
                        for (j = 0; j < nallhits[i]; j++) {
                            fprintf (params->fpout, "%zu %d %.7g\n", qID + i + 1,
                                    allhits[i][j].key + 1, allhits[i][j].val);
                        }
                    }
            fflush (params->fpout);
        }
        params->nSimPairs += da_isum (params->nqrows, nallhits, 1);

        if (params->verbosity > 0 && (qID / params->nqrows) % progressInd == 0)
            da_progress_advance (pct);

    }
    if (params->verbosity > 0) {
        da_progress_finalize (pct);
        printf ("\n");
    }
    params->nDotProducts = params->nCandidates;

    /* free memory */
    da_FreeMatrix ((void ***) &hits, params->nthreads, params->ndrows);
    da_FreeMatrix ((void ***) &cand, params->nthreads, params->ndrows);
    da_FreeMatrix ((void ***) &allhits, params->nqrows, 2 * params->ndrows);
    da_iFreeMatrix (&marker, params->nthreads, params->ndrows);
    da_free ((void **) &nallhits, &nallcands, LTERM);
}




/*************************************************************************/
/*! Computes the neighbors of a set of rows against the documents in
    pmat using OpenMP */
/**************************************************************************/
void da_ompijnnFindNeighbors (params_t * params, da_csr_t * mat, da_csr_t * pmat,
        int qID, int nqrows,
        int dID, int ndrows,
        idx_t * nallhits, idx_t * nallcands,
        da_ivkv_t ** allhits, da_ivkv_t ** atcand,
        da_ivkv_t ** athits, idx_t ** atmarker, char simtype)
{

#pragma omp parallel
    {
        int i, j, k, l, ci, nchits, nhits, ncands, nnbrs, noldhits, tid;
        idx_t *marker;
        da_ivkv_t *cand, *hits, *newhits, *oldhits;

        tid = omp_get_thread_num ();

        marker = atmarker[tid];
        cand = atcand[tid];
        hits = athits[tid];

#pragma omp for private(j, k, l, nhits, nnbrs, newhits, oldhits, noldhits) schedule(dynamic, NITERS)
        for (i = 0; i < nqrows; i++) {
            /* compute the similarity */
            nhits = da_ApssGetSimilarRows (pmat, qID - dID + i, 0 /* include qid==cid case */,
                    mat->rowptr[qID + i + 1] -
                    mat->rowptr[qID + i],
                    mat->rowind + mat->rowptr[qID + i],
                    mat->rowval + mat->rowptr[qID + i],
                    simtype, ndrows, params->eps, hits,
                    marker, cand, &nallcands[i]);

            for (k = 0; k < nhits; k++)
                hits[k].key += dID;

            /* merge with the current best neighbors */
            da_ivkvsortd (nhits, hits);
            nnbrs = da_min (ndrows, nhits);

            newhits = allhits[i];
            oldhits = newhits + nhits;
            noldhits = nallhits[i];
            memmove (oldhits, newhits, sizeof (da_ivkv_t) * noldhits);

            /* the two lists to be merged are (nnbrs, hits) and (noldhits, oldhits) */
            for (l = 0, j = 0, k = 0; j < nnbrs && k < noldhits; l++) {
                if (hits[j].val >= oldhits[k].val) {
                    newhits[l] = hits[j];
                    j++;
                } else {
                    newhits[l] = oldhits[k];
                    k++;
                }
            }
            for (; j < nnbrs; l++, j++)
                newhits[l] = hits[j];
            for (; k < noldhits; l++, k++)
                newhits[l] = oldhits[k];
            nallhits[i] = l;

        }

    }

}



/**
 * Pre-process dataset by pruning, scaling, filtering, shuffling, etc.
 * pre-processing must be done in concert for both docs and docsdb.
 */
void da_ijnnPreProcessData (params_t * params)
{
    /* scale term values */
    if (params->scale) {
        if (params->verbosity > 0)
            printf ("   Scaling input matrices docs and docsdb.\n   If docs is a subset, scaling should be done before split.\n");
        da_csr_Scale (params->docs, DA_SCALE_IDF);
        da_csr_Scale (params->docsdb, DA_SCALE_IDF);
    }

    da_csr_t *docs =  da_csr_Join(params->docs, params->docsdb);
    da_csr_t *tmp = NULL;
    idx_t nrows = params->docs->nrows;
    idx_t nrowsdb = params->docsdb->nrows;
    da_csr_FreeAll(&params->docs, &params->docsdb, LTERM);

    /* prune cols with less than pcminlen or more than pcmaxlen values -- must be done in concert */
    if (params->prunecols) {
        da_errexit("Prunecols not implemented for this mode.");
    }

    /* prune rows with less than prminlen or more than prmaxlen values */
    if (params->prunerows) {
        da_errexit("Prunecols not implemented for this mode.");
    }

    // compact the column space
    da_csr_CompactColumns (docs);

    /* normalize docs rows */
    if (params->norm > 0)
        da_csr_Normalize (docs, DA_ROWCOL, params->norm);

    /* split docs from docsdb */
    params->docs = da_csr_ExtractSubmatrix(docs, 0, nrows);
    params->docsdb = da_csr_ExtractSubmatrix(docs, nrows, nrowsdb);
    da_csr_Free(&docs);

}

