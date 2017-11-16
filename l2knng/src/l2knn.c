/*!
 \file  idxjoin.c
 \brief This file contains idxJoin related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
void refreshIndex(da_iapq_t **knng, da_csr_t *docs, ptr_t *rend, ptr_t *cend,
        val_t *rlen, val_t *clen,
        da_ivpq_t *fsim, idx_t *proc, idx_t nproc, val_t simt);
void finalizeBlock(params_t *params, da_iapq_t **knng, da_ivpq_t *rorder,
        da_csr_t *docs, ptr_t *rend, ptr_t *cend, val_t *rlen, val_t *clen,
        val_t *qrow, val_t *qlen, accum_t *accum, val_t *msims,
        idx_t *cands, idx_t *dcands,
        da_csr_t *neigh, da_ivpq_t *fsim, idx_t *proc, idx_t nproc);
void l2knnFindInitialNeighbors(params_t *params, da_csr_t *docs,
        da_iapq_t **knng, idx_t *cands, idx_t *nnbrs, accum_t *accum, val_t *qrow, da_iakv_t *anbrs);
size_t l2knnEnhanceNeighborhood(params_t *params, da_csr_t *docs,
        da_iapq_t **knng, idx_t *cands, val_t *qrow,
        accum_t *accum, da_iakv_t *anbrs, idx_t *nnbrs);
accum_t sddotp(size_t sz, idx_t *ind, val_t *val, val_t *qrow);


/* dot product between sparse and dense vector */
inline accum_t sddotp(size_t sz, idx_t *ind, val_t *val, val_t *qrow)
{
    size_t i;
    accum_t sim;
    for (sim=0.0, i=0; i < sz; ++i)
        sim += val[i] * qrow[ind[i]];

    return sim;
}


/**
 * Main entry point to L2KNNG.
 */
void l2knnFindNeighbors(params_t *params) {

    ssize_t i, j, k, l, s, b, n, nd, nnz, cnt, dps, isz, rid, cid, ncands;
    idx_t nrows, ncols, progressInd, pct, bid;
    ptr_t *rloc=NULL, *rptr, *rend=NULL, *cptr, *cend=NULL;
    idx_t *rind, *cind, *cands=NULL, *dcands=NULL, *csizes=NULL, *cperm=NULL, *nnbrs=NULL, *proc=NULL;
    val_t maxrval, myval, len, simt, qcms, qms, cms, rct, bnd, srbnd, acc;
    val_t *rval, *cval, *clen=NULL, *rlen=NULL, *msims=NULL, *qrow=NULL, *qlen=NULL;
    da_iakv_t *anbrs=NULL;
    da_iapq_t **knng=NULL;
    da_ivpq_t *rorder=NULL;
    da_ivpq_t *fsim=NULL;
    da_iakv_t *qheap, *cheap;
    da_csr_t *docs=NULL, *sdoc=NULL, *neigh=NULL;
    accum_t *accum=NULL;

    // allocate memory
    docs   = params->docs;
    nrows  = docs->nrows;   // num rows
    ncols  = docs->ncols;   // num cols
    nnz    = docs->rowptr[nrows];
    rval   = docs->rowval;
    rind   = docs->rowind;
    rptr   = docs->rowptr;
    nnbrs  = da_ismalloc(nrows, 0, "l2knnFindNeighbors: nnbrs");  /* number of potential neighbors <= k */
    cands  = da_imalloc(nrows, "l2knnFindNeighbors: candidates"); /* candidates for similarity computation with query row */
    qrow   = da_vsmalloc(ncols, 0, "l2knnFindNeighbors: qrow"); /* dense version of the query vector */
    accum  = da_asmalloc(nrows, -1, "l2knnFindNeighbors: accum"); /* accumulator array */
    csizes = da_ismalloc(ncols, 0.0, "l2knnFindNeighbors: csizes"); /* column sizes */
    cperm  = da_imalloc(ncols, "l2knnFindNeighbors: cperm"); /* permutation vector for columns */
    anbrs  = da_iakvmalloc(da_max(params->k, params->k * params->alpha), "l2knnFindNeighbors: anbrs"); /* initial temporary neighbors */
    rloc   = da_psmalloc(nrows, -1, "l2knnFindNeighbors: rloc"); /* shared locator for all the priority queues */
    proc   = da_imalloc(nrows, "l2knnFindNeighbors: proc");  /* items that have been processed */
    knng   = params->knng = (da_iapq_t**) da_malloc(nrows * sizeof(da_iapq_t*),
            "knng"); /* KNNG result as set of priority queues */
    for (i=0; i < nrows; ++i)
        knng[i] = da_iapqCreateShared(params->k, nrows, rloc);

    /* pre-process data */
    preProcessData(params);

    if(params->verbosity > 0)
        printf("Docs matrix: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, "
            PRNT_PTRTYPE " nnz\n", docs->nrows, docs->ncols, docs->rowptr[docs->nrows]);

    /*** STEP1: Initial KNNG selection based on large values + find rows with less than k neighbors. ***/
    da_startwctimer(params->timer_5); // initial neighbors time
    l2knnFindInitialNeighbors(params, docs, knng, cands, nnbrs, accum, qrow, anbrs);
    da_stopwctimer(params->timer_5); // initial neighbors time
    da_printTimerLong("\tInitial knng build time: ", da_getwctimer(params->timer_5));
    params->timer_3 += params->timer_5;
    da_clearwctimer(params->timer_5);

    /* verify approximate results */
    if (params->vFile)
        verify_knng_results2(knng, nrows, params->vFile, 1);

    /*** STEP2: Enhance initial neighbors with neighbor's neighbors. ***/
    if(params->nenhance){
        for (i=0; i < params->nenhance; ++i) {
            da_startwctimer(params->timer_5); // enhance neighbors time
            cnt = l2knnEnhanceNeighborhood(params, docs, knng, cands, qrow,
                    accum, anbrs, nnbrs);
            printf("# improvements: %zu\n", cnt);
            da_stopwctimer(params->timer_5); // enhance neighbors time
            if (cnt < 0.0001 * nrows * params->k)
                break;
        }
        da_printTimerLong("\tEnhance initial knng time: ", da_getwctimer(params->timer_5));
        params->timer_3 += params->timer_5;
        da_clearwctimer(params->timer_5);

        /* verify approximate results */
        if (params->vFile)
            verify_knng_results2(knng, nrows, params->vFile, 1);
    }

    if(params->mode == MODE_L2KNN_A){  /* only interested in the approximate solution */

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
        if(neigh)
            da_csr_Free(&neigh);
        da_free((void**) &knng, &cperm, &csizes, LTERM);
        da_free((void**) &rloc, &cands, &dcands, &qrow, &qlen, &accum,
                &rlen, &clen, &cend, &msims, &anbrs, &rend, &nnbrs, LTERM);

        return;
    }

    /*** STEP4: Verify found neighbors, adjusting lists as appropriate. ***/
#ifdef L2CV
    qlen = da_vsmalloc(ncols, 0, "l2knnFindNeighbors: qlen"); /* dense hash of suffix l2-norms of x */
#endif
    rorder = da_ivpqCreate(nrows, nrows);
    fsim   = da_ivpqCreate(nrows, nrows);
    msims  = da_vmalloc(nrows, "l2knnFindNeighbors: msims");  /* dense cache of min sim for each row as it gets indexed */
    dcands = da_imalloc(nrows, "l2knnFindNeighbors: dcands"); /* done (pre-computed) candidates for query row */
    rlen   = da_vmalloc(nnz, "l2knnFindNeighbors: rlen");     /* row suffix lengths */
    clen   = da_vmalloc(nnz, "l2knnFindNeighbors: clen");     /* row suffix lengths in inverted index */
    rend   = da_pmalloc(nrows, "l2knnFindNeighbors: rend");   /* end pointers for the partial dynamic forward index */
    cend   = da_pmalloc(ncols, "l2knnFindNeighbors: cend");   /* end pointers for the partial dynamic inverted index */
    cptr   = docs->colptr;  /* reusing space from col structure for inverted index */
    cind   = docs->colind;
    cval   = docs->colval;

    da_startwctimer(params->timer_5);
    if (params->verbosity > 0) {
        printf("Verifying ");
        fflush(stdout);
    }
    /* reorder columns, reset column pointers & copy them to column end pointers too */
    l2knnReorderCols(docs, csizes, cperm);
    for(n=0, i=0; i < ncols; ++i){
        cptr[i] = cend[i] = n;
        n += csizes[i];
    }
    da_free((void**)&cperm, &csizes, LTERM);
    params->cperm = NULL;

    /* set up index pointers for the forward index */
    da_pcopy(nrows, rptr, rend);

    /* create inverted index of current found neighbors */
    neigh = knng2csr(knng, nrows);
    da_csr_CreateIndex(neigh, DA_COL);

    /* find order of row processing */
    for (i=0; i < nrows; ++i) {
        if(nnbrs[i] >= params->k)
            da_ivpqInsert(rorder, i, da_iapqSeeTopVal(knng[i]));
    }
    isz = ceil(rorder->nnodes/(1.0*params->nbl));
    if (params->verbosity > 0) {
        printf("%zu neighborhoods (in blocks of %zu)... ", rorder->nnodes, isz);
        fflush(stdout);
    }

    /* execute search */
    n = 0;
    da_progress_init(pct, progressInd, rorder->nnodes);
    while(rorder->nnodes > 0){
        rid = rorder->heap[0].key;
        rct = da_ivpqSeeTopVal(fsim);     /* threshold for accepting new candidates */
        qheap = knng[rid]->heap;
        qms = qheap[0].val;               /* query vector minimum neighborhood similarity */

        da_iapqInitLocator(knng[rid]);

        #ifdef EXTRATIMES
        da_startwctimer(params->timer_1);  /* candidate generation time */
        #endif

        /* mark items that we already have answers for:
         *      first current neighborhood, than reverse neighborhood */
        nd = 0;
        for(i=0; i < knng[rid]->nnodes; ++i){
            cid = qheap[i].key;
            accum[cid] = -2;
            dcands[nd++] = cid;
        }
        for(i=neigh->colptr[rid]; i < neigh->colptr[rid+1]; i++){
            cid = neigh->colind[i];
            accum[cid] = -2;
            dcands[nd++] = cid;
        }

        /*** find matches for query row ***/
        for (bnd=1.0, srbnd=1.0, ncands=0, i=rptr[rid]; i < rptr[rid+1]; ++i) {
            myval = rval[i];
            bnd -= myval * myval;

            /* cache the query row */
            k = rind[i];
            qrow[k] = myval;
#ifdef L2CV
            qlen[k] = bnd > 0 ? sqrt(bnd) : 0.0;
#endif

            for (j=cptr[k]; j < cend[k]; ++j) {
                cid = cind[j];  /* potential candidate */
                acc = accum[cid];
                if (acc > 0.0) {
                    // accumulate
                    acc = accum[cid] += myval * cval[j];
#ifdef L2CG
                    // check l2-norm bound
                    if (acc + srbnd * clen[j] < da_min(qms, msims[cid])) {
                        accum[cid] = -2;
                        #ifdef EXTRACOUNTS
                        params->nPruneLength++;
                        #endif
                    }
#endif
                } else if (srbnd >= rct && acc > -2) {
                    // accumulate
                    if (acc == -1) {
                        acc = accum[cid] = myval * cval[j];
                        cands[ncands++] = cid;
                    } else
                        acc = accum[cid] += myval * cval[j];
#ifdef L2CG
                    // check l2-norm bound
                    if (acc + srbnd * clen[j] < da_min(qms, msims[cid])) {
                        accum[cid] = -2;
                        #ifdef EXTRACOUNTS
                        params->nPruneLength++;
                        #endif
                    }
#endif
                }
            }


            srbnd = bnd > 0 ? sqrt(bnd) : 0.0;
        }


        #ifdef EXTRATIMES
        da_stopwctimer(params->timer_1);  /* candidate generation time */
        da_startwctimer(params->timer_2); /* candidate verification time */
        #endif

        /*** process candidates ***/
        params->nCandidates += ncands;
        for (dps=0, i=0; i < ncands; i++) {
            cid = cands[i];
            if (accum[cid] <= 0.0) {
                accum[cid] = -1;
                continue;
            }

            cheap = knng[cid]->heap;
            acc = accum[cid];
            accum[cid] = -1;
            j = rend[cid];
#if defined(PSCV) || defined(L2CV)
            qcms = da_min(qms, msims[cid]);   /* min sim of either the query or candidate neighborhoods */
#endif

#ifdef PSCV
            if (j < rptr[cid+1] && acc + rlen[j] < qcms && qrow[rind[j]] == 0) {
                #ifdef EXTRACOUNTS
                params->nPrunePscore++;
                #endif
                continue;
            }
#endif

            /* compute the remaining portion of the dot product */
            for(; j < rptr[cid+1]; ++j){
                if(qrow[rind[j]] > 0){
                    acc += rval[j] * qrow[rind[j]];
                    /* check L2-norm bound */
#ifdef L2CV
                    if(acc + rlen[j] * qlen[rind[j]] < qcms){
                        #ifdef EXTRACOUNTS
                        params->nPruneLength2++;
                        #endif
                        goto nexcid;
                    }
#endif
                }
            }

            dps++;

            if(acc > qms){
                da_iapqInsert(knng[rid], cid, acc);
                qms = da_iapqSeeTopVal(knng[rid]);
            }
            if (acc > cheap[0].val) {
                cms = cheap[0].val;
                da_iapqInsertHeap(knng[cid], rid, acc);
                if(cheap[0].val > cms){
                    msims[cid] = cheap[0].val;
                    da_ivpqUpdate(fsim, cid, msims[cid]);
                }
            }

            nexcid:
                continue;
        }

        params->nDotProducts += dps;

        /* clear the accumulation array for the pre-computed candidates */
        for(i=0; i < nd; ++i)
            accum[dcands[i]] = -1;

        #ifdef EXTRATIMES
        da_stopwctimer(params->timer_2);  /* candidate verification time */
        #endif

        memcheckonly:

        /* reset the heap locator for the current row */
        da_iapqResetLocator(knng[rid]);

        #ifdef EXTRATIMES
        da_startwctimer(params->timer_4); /* indexing time */
        #endif

        /* set the row similarity at time of indexing */
        simt = rorder->heap[0].val;       /* query vector indexing similarity threshold */
        da_ivpqGetTop(rorder);
        da_ivpqInsert(fsim, rid, simt);
        msims[rid] = qms;

        /* un-cache the query row & add forward index data */
        for(bnd=1.0, rend[rid]=rptr[rid+1], j=rptr[rid]; j < rptr[rid+1]; ++j){
            qrow[rind[j]] = 0;

            /** index part of the row **/
            bnd -= rval[j] * rval[j];
            srbnd = bnd > 0 ? sqrt(bnd) : 0.0;
            b = cend[rind[j]]++;
            cval[b] = rval[j];
            cind[b] = rid;
            clen[b] = srbnd;

            /** check where to stop indexing **/
            if(srbnd < simt){
                rend[rid] = ++j;
                break;
            }
        }
        for(; j < rptr[rid+1]; ++j){
            qrow[rind[j]] = 0;
            bnd -= rval[j] * rval[j];
            rlen[j] = bnd > 0 ? sqrt(bnd) : 0.0;
        }

        #ifdef EXTRATIMES
        da_stopwctimer(params->timer_4); /* indexing time */
        #endif

        if (params->verbosity > 0 && rorder->nnodes % progressInd == 0){
            da_progress_advance(pct);
            printf(" %.2f ", simt);
        }

        proc[n++] = rid;  /* rid finished processing */

        if(n == isz && rorder->nnodes >= n){
            printf("f ");
            fflush(stdout);
            finalizeBlock(params, knng, rorder, docs, rend, cend, rlen, clen, qrow, qlen, accum, msims,
                    cands, dcands, neigh, fsim, proc, n);
            n=0;
        }

    }

    if (params->verbosity > 0)
        da_progress_finalize(pct);
    printf("\n");
    da_stopwctimer(params->timer_5);
    da_printTimerLong("\tVerification time: ", da_getwctimer(params->timer_5));
    params->timer_3 += params->timer_5;

    /* verify results */
    if (params->vFile)
        verify_knng_results2(knng, nrows, params->vFile, params->verbosity);

    /* report number of similarity pairs found */
    for(i=0; i < nrows; ++i)
        params->nSimPairs += knng[i]->nnodes;

    /* save neighbor graph */
    da_save_knng(params, knng, nrows, DA_FMT_CLUTO);

    /* find final inverted index size */
    for(n=0, i=0; i < ncols; ++i)
        n += cend[i] - cptr[i];
    params->indexSize += n;

    /* free memory */
    for (i=0; i < nrows; ++i) {
        knng[i]->locator = NULL;
        da_iapqDestroy(knng[i]);
    }
    if(rorder)
        da_ivpqDestroy(rorder);
    if(fsim)
        da_ivpqDestroy(fsim);
    if(neigh)
        da_csr_Free(&neigh);
    da_free((void**) &knng, LTERM);
    da_free((void**) &rloc, &cands, &dcands, &qrow, &qlen, &accum, &proc,
            &rlen, &clen, &cend, &msims, &anbrs, &rend, &nnbrs, LTERM);

}


/**
 * Finalize neighborhoods for partially executed rows in proc and reset the inverted index
 */
void finalizeBlock(params_t *params, da_iapq_t **knng, da_ivpq_t *rorder,
        da_csr_t *docs, ptr_t *rend, ptr_t *cend, val_t *rlen, val_t *clen,
        val_t *qrow, val_t *qlen, accum_t *accum, val_t *msims,
        idx_t *cands, idx_t *dcands,
        da_csr_t *neigh, da_ivpq_t *fsim, idx_t *proc, idx_t nproc)
{

    size_t i, j, k, ri, rid, cid, nd, ncands, dps;
    val_t bnd, srbnd, qms, cms, qcms, rct, myval, acc;
    ptr_t b, e;
    ptr_t *rptr, *cptr;
    idx_t *rind, *cind;
    val_t *rval, *cval;
    da_iakv_t *qheap, *cheap;

    rptr   = docs->rowptr;
    rind   = docs->rowind;
    rval   = docs->rowval;
    cptr   = docs->colptr;  /* reusing space from col structure for inverted index */
    cind   = docs->colind;
    cval   = docs->colval;

    /** add the current inverted index size **/
    for(k=0, i=0; i < docs->ncols; ++i)
        k += cend[i] - cptr[i];
    params->indexSize += k;

    for(ri=0; ri < rorder->nnodes; ++ri){
        rid = rorder->heap[ri].key;
        rct = da_ivpqSeeTopVal(fsim);     /* threshold for accepting new candidates */
        qheap = knng[rid]->heap;
        qms = qheap[0].val;               /* query vector minimum neighborhood similarity */

        da_iapqInitLocator(knng[rid]);

        #ifdef EXTRATIMES
        da_startwctimer(params->timer_1);  /* candidate generation time */
        #endif

        /* mark items that we already have answers for:
         *      first current neighborhood, than reverse neighborhood */
        nd = 0;
        for(i=0; i < knng[rid]->nnodes; ++i){
            cid = qheap[i].key;
            accum[cid] = -2;
            dcands[nd++] = cid;
        }
        for(i=neigh->colptr[rid]; i < neigh->colptr[rid+1]; i++){
            cid = neigh->colind[i];
            accum[cid] = -2;
            dcands[nd++] = cid;
        }

        /*** find matches for query row ***/
        for (bnd=1.0, srbnd=1.0, ncands=0, i=rptr[rid]; i < rptr[rid+1]; ++i) {
            myval = rval[i];
            bnd -= myval * myval;

            /* cache the query row */
            k = rind[i];
            qrow[k] = myval;
#ifdef L2CV
            qlen[k] = bnd > 0 ? sqrt(bnd) : 0.0;
#endif

            for (j=cptr[k]; j < cend[k]; ++j) {
                cid = cind[j];  /* potential candidate */
                acc = accum[cid];
                if (acc > 0.0) {
                    // accumulate
                    acc = accum[cid] += myval * cval[j];
#ifdef L2CG
                    // check l2-norm bound
                    if (acc + srbnd * clen[j] < da_min(qms, msims[cid])) {
                        accum[cid] = -2;
                        #ifdef EXTRACOUNTS
                        params->nPruneLength++;
                        #endif
                    }
#endif
                } else if (srbnd >= rct && acc > -2) {
                    // accumulate
                    if (acc == -1) {
                        acc = accum[cid] = myval * cval[j];
                        cands[ncands++] = cid;
                    } else
                        acc = accum[cid] += myval * cval[j];
#ifdef L2CG
                    // check l2-norm bound
                    if (acc + srbnd * clen[j] < da_min(qms, msims[cid])) {
                        accum[cid] = -2;
                        #ifdef EXTRACOUNTS
                        params->nPruneLength++;
                        #endif
                    }
#endif
                }
            }


            srbnd = bnd > 0 ? sqrt(bnd) : 0.0;
        }

        #ifdef EXTRATIMES
        da_stopwctimer(params->timer_1);  /* candidate generation time */
        da_startwctimer(params->timer_2); /* candidate verification time */
        #endif

        /*** process candidates ***/
        params->nCandidates += ncands;
        for (dps=0, i=0; i < ncands; i++) {
            cid = cands[i];
            if (accum[cid] <= 0.0) {
                accum[cid] = -1;
                continue;
            }

            cheap = knng[cid]->heap;
            acc = accum[cid];
            accum[cid] = -1;
            j = rend[cid];
#if defined(PSCV) || defined(L2CV)
            qcms = da_min(qms, msims[cid]);   /* min sim of either the query or candidate neighborhoods */
#endif

#ifdef PSCV
            if (j < rptr[cid+1] && acc + rlen[j] < qcms && qrow[rind[j]] == 0) {
                #ifdef EXTRACOUNTS
                params->nPrunePscore++;
                #endif
                continue;
            }
#endif

            /* compute the remaining portion of the dot product */
            for(; j < rptr[cid+1]; ++j){
                if(qrow[rind[j]] > 0){
                    acc += rval[j] * qrow[rind[j]];
                    /* check L2-norm bound */
#ifdef L2CV
                    if(acc + rlen[j] * qlen[rind[j]] < qcms){
                        #ifdef EXTRACOUNTS
                        params->nPruneLength2++;
                        #endif
                        goto nexcid;
                    }
#endif
                }
            }

            dps++;

            if(acc > qms){
                da_iapqInsert(knng[rid], cid, acc);
                qms = da_iapqSeeTopVal(knng[rid]);
            }
            if (acc > cheap[0].val) {
                cms = cheap[0].val;
                da_iapqInsertHeap(knng[cid], rid, acc);
                if(cheap[0].val > cms){
                    msims[cid] = cheap[0].val;
                    da_ivpqUpdate(fsim, cid, msims[cid]);
                }
            }

            nexcid:
                continue;
        }

        params->nDotProducts += dps;

        /* reset the query row */
        for(i=rptr[rid]; i < rptr[rid+1]; ++i)
            qrow[rind[i]] = 0.0;


        /* clear the accumulation array for the pre-computed candidates */
        for(i=0; i < nd; ++i)
            accum[dcands[i]] = -1;

        #ifdef EXTRATIMES
        da_stopwctimer(params->timer_2);  /* candidate verification time */
        #endif

        /* reset the heap locator for the current row */
        da_iapqResetLocator(knng[rid]);

    }

    /* remove done rows from fsim */
    da_ivpqReset(fsim);

    /* reset the inverted index */
    memcpy(cend, cptr, docs->ncols * sizeof(ptr_t));

    /* set new updated order and sim thresholds for the remaining rows */
    k = rorder->nnodes;
    for(i=0; i < k; ++i)
        cands[i] = rorder->heap[i].key;
    da_ivpqReset(rorder);
    for(i=0; i < k; ++i)
        da_ivpqInsert(rorder, cands[i], da_iapqSeeTopVal(knng[cands[i]]));

}



/**
 * Refresh the index given updated sim values
 */
void refreshIndex(da_iapq_t **knng, da_csr_t *docs, ptr_t *rend, ptr_t *cend,
        val_t *rlen, val_t *clen,
        da_ivpq_t *fsim, idx_t *proc, idx_t nproc, val_t simt)
{
    size_t i, j, rid;
    val_t bnd, srbnd, qms;
    ptr_t b, e;
    ptr_t *rptr, *cptr;
    idx_t *rind, *cind;
    val_t *rval, *cval;

    rptr   = docs->rowptr;
    rind   = docs->rowind;
    rval   = docs->rowval;
    cptr   = docs->colptr;  /* reusing space from col structure for inverted index */
    cind   = docs->colind;
    cval   = docs->colval;

    /* reset index pointers */
    memcpy(cend, cptr, docs->ncols * sizeof(ptr_t));

    for(i=0; i < nproc; ++i){
        rid = proc[i];
        qms = da_iapqSeeTopVal(knng[rid]);
        qms = da_min(qms, simt);
        da_ivpqUpdate(fsim, rid, qms);
        e = rend[rid];
        for(bnd=1.0, rend[rid]=rptr[rid+1], j=rptr[rid]; j < rend[rid]; ++j){
            /** index part of the row **/
            bnd -= rval[j] * rval[j];
            srbnd = bnd > 0 ? sqrt(bnd) : 0.0;
            b = cend[rind[j]]++;
            cval[b] = rval[j];
            cind[b] = rid;
            clen[b] = srbnd;

            /** check where to stop indexing **/
            if(srbnd < qms){
                rend[rid] = ++j;
                break;
            }
        }
        for(; j < e; ++j){
            bnd -= rval[j] * rval[j];
            rlen[j] = bnd > 0 ? sqrt(bnd) : 0.0;
        }
    }
}

size_t l2knnEnhanceNeighborhood(params_t *params, da_csr_t *docs,
        da_iapq_t **knng, idx_t *cands, val_t *qrow,
        accum_t *accum, da_iakv_t *anbrs, idx_t *nnbrs) {
    size_t i, j, k, l, nn, ncands, nadded = 0;
    idx_t nrows, cid, progressInd, pct;
    val_t rv, simv, minsim;
    da_csr_t *neigh;
    da_iapq_t *lnn;
    ptr_t *rowptr, *colptr, *ptr;
    idx_t *rowind, *colind, *ind;
    val_t *rowval, *colval, *val;
    char *tag;

    nrows = docs->nrows; /* num rows */
    nn    = params->k * params->alpha; /* number of neighbors' neighbors to retrieve */
    neigh = knng2csr(knng, nrows);
    da_csr_CreateIndex(neigh, DA_COL);
    da_csr_SortValues(neigh, DA_COL, DA_SORT_D);
    da_csr_SortValues(neigh, DA_ROW, DA_SORT_D);

    ptr = docs->rowptr;
    ind = docs->rowind;
    val = docs->rowval;

    rowptr = neigh->rowptr;
    rowind = neigh->rowind;
    rowval = neigh->rowval;
    tag = da_csmalloc(nrows, 0, "l2knnEnhanceNeighborhood: tag");

    if (params->verbosity > 0) {
        printf("Enhance neighborhood... ");
    }
    da_progress_init(pct, progressInd, nrows);

    colptr = neigh->colptr;
    colind = neigh->colind;
    colval = neigh->colval;

    for (i=0; i < nrows; ++i) {
        if (rowptr[i+1] == rowptr[i] || nnbrs[i] < params->k)
            continue; /*  no neighbors or already found all the neighbors */
        lnn = knng[i];
        da_iapqInitLocator(lnn);
        for (j=ptr[i]; j < ptr[i+1]; ++j) /* cache query row */
            qrow[ind[j]] = val[j];

        /* tag current neighbors */
        for(j=rowptr[i]; j < rowptr[i+1]; ++j)
            tag[rowind[j]] = 1;

        /* find neighbor's neighbors */
        for (ncands=0, j=rowptr[i]; j < rowptr[i+1]; ++j) {
            cid = rowind[j];

            for (k=rowptr[cid]; k < rowptr[cid+1]; ++k) { /* neighbors's neighbors */
#ifdef PRUNENGBS
                if(rowval[k] < rowval[j])
                    break;
#endif
                if(ncands >= nn)
                    break;
                if (rowind[k] != i && tag[rowind[k]] == 0) {
                    cands[ncands++] = rowind[k];
                    tag[rowind[k]] = 1;
                }
            }
        }

        for (k=colptr[i]; k < colptr[i]; ++k) { /* reverse neighbors */
#ifdef PRUNENGBS
            if(colval[k] < rowval[rowptr[i+1]-1])
                break;
#endif
            if(ncands >= nn)
                break;
            if (tag[colind[k]] == 0) {
                cands[ncands++] = colind[k];
                tag[colind[k]] = 1;
            }
        }

        /* untag current neighbors */
        for(j=rowptr[i]; j < rowptr[i+1]; ++j)
            tag[rowind[j]] = 0;


        minsim = da_iapqSeeTopVal(lnn);
        for (j=0; j < ncands; ++j) {
            cid = cands[j];
            tag[cid] = 0;
#ifdef EXTRACOUNTS
            params->nDotProducts2++;
#endif
            simv = sddotp(ptr[cid+1] - ptr[cid], ind + ptr[cid],
                    val + ptr[cid], qrow);
            if (simv > minsim) {
                da_iapqInsert(lnn, cid, simv);
                minsim = da_iapqSeeTopVal(lnn);
                nadded++;
            }
            if (simv > da_iapqSeeTopVal(knng[cid])) {
                if(da_iapqInsertHeap(knng[cid], i, simv))
                    nadded++;
            }
        }

        da_iapqResetLocator(lnn);
        for (j=ptr[i]; j < ptr[i+1]; ++j)
            qrow[ind[j]] = 0; /* reset query row */

        if (params->verbosity > 0 && i % progressInd == 0)
            da_progress_advance(pct);
    }
    if (params->verbosity > 0)
        da_progress_finalize(pct);
    printf("\n");

    da_csr_Free(&neigh);
    da_free((void**)&tag, LTERM);
    return nadded;
}

/**
 * Find initial neighbors that are likely to be true neighbors
 */
void l2knnFindInitialNeighbors(params_t *params, da_csr_t *docs,
        da_iapq_t **knng, idx_t *cands, idx_t *nnbrs, accum_t *accum,
        val_t *qrow, da_iakv_t *anbrs) {
    ssize_t i, j, k, nn, c;
    size_t nns, nns2, lj;
    idx_t r, c1, c2, cid, nrows, ncols, progressInd, pct;
    ptr_t cs1, ce1, cs2, ce2, *ptr;
    val_t cv1, cv2, simv, v;
    ptr_t *rowptr, *colptr;
    idx_t *rowind, *colind;
    val_t *rowval, *colval;
    da_ivpq_t *mq = NULL;

    nrows = docs->nrows;   // num rows
    ncols = docs->ncols;   // num cols
    nn = params->k * params->alpha; /* number of neighbor's neighbors to retrieve */

    if(params->alpha < 1){  /* randomly choose the initial k neighbors for each row & identify potential neighborhood size */
        if (!docs->colptr)
            da_csr_CreateIndex(docs, DA_COL);
        rowptr = docs->rowptr; /* index in rowind/rowval where each row starts */
        rowind = docs->rowind; /* column ids in each row */
        rowval = docs->rowval; /* values in each row */
        colptr = docs->colptr; /* index in rowind/rowval where each row starts */
        colind = docs->colind; /* column ids in each col */
        colval = docs->colval; /* values in each col */
        nn     = params->k;

        if (params->verbosity > 0) {
            printf("nn: %zu\n", nn);
            printf("Find initial neighbors... ");
        }

        da_progress_init(pct, progressInd, nrows);
        for (i=0; i < nrows; ++i) {
            if (rowptr[i+1] == rowptr[i])
                continue; /* empty row, no neighbors */
            /* find out potential number of neighbors */
            for(nns=0, j=rowptr[i]; j < rowptr[i+1] && nns < nn; ++j){
                for(k=colptr[rowind[j]]; k < colptr[rowind[j]+1] && nns < nn; ++k){
                    if(colind[k] != i && accum[colind[k]] < 0){
                        cands[nns++] = colind[k];
                        accum[colind[k]] = 1;
                    }
                }
            }
            nnbrs[i] = nns;
            for(j=0; j < nns; ++j)
                accum[cands[j]] = -1;
            if(nns == 0)
                continue;
            /* randomly choose initial neighborhood */
            da_iapqInitLocator(knng[i]);
            for(j=rowptr[i]; j < rowptr[i+1]; ++j)
                qrow[rowind[j]] = rowval[j];
            if(nns >= nn){
                nns2 = 0;
                do {
                    j = rowptr[i] + da_irandInRange(rowptr[i+1]-rowptr[i]);
                    k = colptr[rowind[j]] + da_irandInRange(colptr[rowind[j]+1]-colptr[rowind[j]]);
                    if(colind[k] == i)
                        continue;
                    cid = colind[k];
                    simv = sddotp(rowptr[cid+1] - rowptr[cid], rowind + rowptr[cid],
                            rowval + rowptr[cid], qrow);
                    /* insert into graph */
                    da_iapqInsert(knng[i], cid, simv);
                    if (knng[cid]->nnodes < params->k || simv > da_iapqSeeTopVal(knng[cid]))
                        da_iapqInsertHeap(knng[cid], i, simv);
                    nns2++;
                } while(nns2 < nns);
            } else {  /* row has less than k potential neighbors - compute their similarities */
                for(j=0; j < nns; ++j){
                    cid = cands[j];
                    simv = sddotp(rowptr[cid+1] - rowptr[cid], rowind + rowptr[cid],
                            rowval + rowptr[cid], qrow);
                    /* insert into graph */
                    da_iapqInsert(knng[i], cid, simv);
                    if (knng[cid]->nnodes < params->k || simv > da_iapqSeeTopVal(knng[cid])){
                        da_iapqInsertHeap(knng[cid], i, simv);
                    }
                }
            }
            for(j=rowptr[i]; j < rowptr[i+1]; ++j)
                qrow[rowind[j]] = 0.0;
            da_iapqResetLocator(knng[i]);
            #ifdef EXTRACOUNTS
            params->nDotProducts1 += nns;
            #endif

            if (params->verbosity > 0 && i % progressInd == 0)
                da_progress_advance(pct);
        }
        if (params->verbosity > 0)
            da_progress_finalize(pct);
        printf("\n");

        return;
    }

    /* build inverted index and sort both rows and columns in decreasing value oder */
    if (!docs->colptr)
        da_csr_CreateIndex(docs, DA_COL);
    da_csr_SortValues(docs, DA_ROW, DA_SORT_D);
    da_csr_SortValues(docs, DA_COL, DA_SORT_D);
    rowptr = docs->rowptr; /* index in rowind/rowval where each row starts */
    rowind = docs->rowind; /* column ids in each row */
    rowval = docs->rowval; /* values in each row */
    colptr = docs->colptr; /* index in rowind/rowval where each row starts */
    colind = docs->colind; /* column ids in each col */
    colval = docs->colval; /* values in each col */

    if (params->verbosity > 0) {
        printf("nn: %zu\n", nn);
        printf("Find initial neighbors... ");
    }

    da_progress_init(pct, progressInd, nrows);
    if (params->pmerge == PMERGE_ALL){

        mq = da_ivmqCreate(nn, ncols);
        ptr = da_pmalloc(ncols, "l2knnFindInitialNeighbors: ptr");

        for (i=0; i < nrows; ++i) {
            if (rowptr[i+1] == rowptr[i])
                continue; /* empty row, no neighbors */

            for(j=rowptr[i]; j < rowptr[i+1]; ++j){
                qrow[rowind[j]] = rowval[j];
                k = colptr[rowind[j]];
                da_ivmqInsert(mq, rowind[j], rowval[j]*colval[k]);
                ptr[rowind[j]] = k;
            }
            nns = 0;
            while(nns < nn){
                v = da_ivmqSeeTopVal(mq);  /* max value in heap */
                c = da_ivmqGetTop(mq);     /* column value came from */
                if(c < 0)
                    break;                 /* out of items */
                r = colind[ptr[c]];         /* row id for that associated value */
                if(r == i)
                    continue;
                if(accum[r] < 0){
                    anbrs[nns++].key = r;
                    accum[r] = v;
                } else
                    accum[r] += v;

                if(++ptr[c] < colptr[c+1])
                    da_ivmqInsert(mq, c, qrow[c]*colval[ptr[c]]);

            }

            if (nns == 0) {
                for (k = rowptr[i]; k < rowptr[i+1]; ++k)
                    qrow[rowind[k]] = 0;
                continue;
            }
            nnbrs[i] = nns;

            for (k=0; k < nns; ++k) {
                anbrs[k].val = accum[anbrs[k].key];
                accum[anbrs[k].key] = -1;
            }
            da_iakvsortd(nns, anbrs);

            /* compute similarities of initially identified neighbors */
            da_iapqInitLocator(knng[i]);
            for (k=0; k < nns; ++k) {
                cid = anbrs[k].key;
                simv = sddotp(rowptr[cid+1] - rowptr[cid], rowind + rowptr[cid],
                        rowval + rowptr[cid], qrow);
                da_iapqInsert(knng[i], cid, simv);
                if (knng[cid]->nnodes < params->k || simv > da_iapqSeeTopVal(knng[cid]))
                    da_iapqInsertHeap(knng[cid], i, simv);
            }
            da_iapqResetLocator(knng[i]);
#ifdef EXTRACOUNTS
            params->nDotProducts1 += nns;
#endif
            for (k=rowptr[i]; k < rowptr[i+1]; ++k)
                qrow[rowind[k]] = 0;

            if (params->verbosity > 0 && i % progressInd == 0)
                da_progress_advance(pct);
        }
        if (params->verbosity > 0)
            da_progress_finalize(pct);
        printf("\n");

        if(mq)
            da_ivmqDestroy(mq);
        da_free((void**)&ptr, &accum, LTERM);

        return;
    }

    for (i=0; i < nrows; ++i) {
        if (rowptr[i + 1] == rowptr[i])
            continue; /* empty row, no neighbors */

        j = rowptr[i];
        lj = rowptr[i+1];
        c1 = rowind[j];
        cv1 = rowval[j];
        qrow[c1] = cv1;
        cs1 = colptr[c1];
        ce1 = colptr[c1+1];
        cv2 = cs2 = ce2 = 0;
        if (j+1 < lj) {
            j++;
            c2 = rowind[j];
            cv2 = rowval[j];
            qrow[c2] = cv2;
            cs2 = colptr[c2];
            ce2 = colptr[c2+1];
        } else
            c2 = -1;
        nns = 0;
        while (nns < nn) {
            if (c2 < 0) { /* last col has been reached. */
                for (; cs1 < ce1 && nns < nn; ++cs1) {
                    if (colind[cs1] == i)
                        continue;
                    if (accum[colind[cs1]] < 0) {
                        anbrs[nns++].key = colind[cs1];
                        accum[colind[cs1]] = colval[cs1] * cv1;
                    } else
                        accum[colind[cs1]] += colval[cs1] * cv1;
                }
                break;
            } else if (c1 < 0) { /* last col has been reached. */
                for (; cs2 < ce2 && nns < nn; ++cs2) {
                    if (colind[cs2] == i)
                        continue;
                    if (accum[colind[cs2]] < 0) {
                        anbrs[nns++].key = colind[cs2];
                        accum[colind[cs2]] = colval[cs2] * cv2;
                    } else
                        accum[colind[cs2]] += colval[cs2] * cv2;
                }
                break;
            }

            if (colval[cs1] * cv1 > colval[cs2] * cv2) {
                if (colind[cs1] != i) {
                    if (accum[colind[cs1]] < 0) {
                        anbrs[nns++].key = colind[cs1];
                        accum[colind[cs1]] = colval[cs1] * cv1;
                    } else
                        accum[colind[cs1]] += colval[cs1] * cv1;

                }
                cs1++;

            } else {
                if (colind[cs2] != i) {
                    if (accum[colind[cs2]] < 0) {
                        anbrs[nns++].key = colind[cs2];
                        accum[colind[cs2]] = colval[cs2] * cv2;
                    } else
                        accum[colind[cs2]] += colval[cs2] * cv2;
                }
                cs2++;
            }

            if (cs1 == ce1) { /* end of col 1 - switch to another col, if exists */
                if (j+1 < lj) {
                    j++;
                    c1 = rowind[j];
                    cv1 = rowval[j];
                    qrow[c1] = cv1;
                    cs1 = colptr[c1];
                    ce1 = colptr[c1 + 1];
                } else
                    c1 = -1;
            } else if (cs2 == ce2) { /* end of col 2 - switch to another col, if exists */
                if (j+1 < lj) {
                    j++;
                    c2 = rowind[j];
                    cv2 = rowval[j];
                    qrow[c2] = cv2;
                    cs2 = colptr[c2];
                    ce2 = colptr[c2 + 1];
                } else
                    c2 = -1;
            }
        }
        if (nns == 0) {
            for (k = rowptr[i]; k < j; ++k)
                qrow[rowind[k]] = 0;
            continue;
        }
        nnbrs[i] = nns;

        for (; j < lj; ++j)
            qrow[rowind[j]] = rowval[j];

        for (k=0; k < nns; ++k) {
            anbrs[k].val = accum[anbrs[k].key];
            accum[anbrs[k].key] = -1;
        }
        da_iakvsortd(nns, anbrs);

        /* compute similarities of initially identified neighbors */
        da_iapqInitLocator(knng[i]);
        for (k=0; k < nns; ++k) {
            cid = anbrs[k].key;
            simv = sddotp(rowptr[cid+1] - rowptr[cid], rowind + rowptr[cid],
                    rowval + rowptr[cid], qrow);
            da_iapqInsert(knng[i], cid, simv);
            if (knng[cid]->nnodes < params->k || simv > da_iapqSeeTopVal(knng[cid]))
                da_iapqInsertHeap(knng[cid], i, simv);
        }
        da_iapqResetLocator(knng[i]);
#ifdef EXTRACOUNTS
        params->nDotProducts1 += nns;
#endif
        for (k=rowptr[i]; k < rowptr[i+1]; ++k)
            qrow[rowind[k]] = 0;

        if (params->verbosity > 0 && i % progressInd == 0)
            da_progress_advance(pct);
    }
    if (params->verbosity > 0)
        da_progress_finalize(pct);
    printf("\n");

}

/**
 * Reorder the document matrix in decreasing number of column nnzs.
 * Also compute and store necessary stats about the data.
 * Store initial order in cperm.
 */
void l2knnReorderCols(da_csr_t *docs, idx_t *csizes, idx_t *cperm) {
    ssize_t i, j, k, nnz;
    idx_t nrows, ncols;
    da_iikv_t *orderColNnz = NULL;
    da_ivkv_t *orderColVal = NULL;
    char pv = 0;
    ptr_t *rowptr;
    idx_t *rowind, *tsizes = NULL;
    val_t v, *rowval, *cval = NULL;

    nrows = docs->nrows;
    ncols = docs->ncols;
    rowptr = docs->rowptr;
    rowind = docs->rowind;
    nnz = rowptr[nrows];
#if COLSORTBY == COLSORTBY_MAX || COLSORTBY == COLSORTBY_AVG
    rowval = docs->rowval;
    cval = da_vsmalloc(ncols, 0.0, "l2knnReorderCols: cval");
#endif

    //record col nnzs
    da_iset(ncols, 0, csizes);
    for (i = 0; i < nnz; i++) {
#if COLSORTBY == COLSORTBY_MAX
        j = rowind[i];
        if(rowval[j] > cval[j])
        cval[j] = rowval[j];
        csizes[j]++;
#elif COLSORTBY == COLSORTBY_AVG
        j = rowind[i];
        cval[j] += rowval[j];
        csizes[j]++;
#else
        csizes[rowind[i]]++;
#endif
    }

#if COLSORTBY == COLSORTBY_AVG
    for(i=0; i < ncols; ++i)
        if(csizes[i] > 0)
            cval[i] /= csizes[i];
#endif

#if COLSORTBY == COLSORTBY_MAX || COLSORTBY == COLSORTBY_AVG
    orderColVal = da_ivkvmalloc(ncols, "l2knnReorderCols: orderColNnz");
    //get new column order
    for (j=0; j < ncols; j++) {
        orderColVal[j].key = j;
        orderColVal[j].val = cval[j];
    }
    da_ivkvsorti(ncols, orderColVal); // sort columns

    for (j=0; j < ncols; j++) {
        cperm[orderColVal[j].key] = j;  // store column permutation
        csizes[j] = orderColVal[j].val;// new column size after permuting
    }
#elif COLSORTBY == COLSORTBY_SIZE
    orderColNnz = da_iikvmalloc(ncols, "l2knnReorderCols: orderColNnz");
    //get new column order
    for (j=0; j < ncols; j++) {
        orderColNnz[j].key = j;
        orderColNnz[j].val = csizes[j];
    }
    da_iikvsorti(ncols, orderColNnz); // sort columns

    for (j=0; j < ncols; j++) {
        cperm[orderColNnz[j].key] = j;  // store column permutation
        csizes[j] = orderColNnz[j].val; // new column size after permuting
    }
#elif COLSORTBY == COLSORTBY_RAND
    da_irandArrayPermute(ncols, cperm, 10*ncols, 1);
    tsizes = da_imalloc(ncols, "l2knnReorderCols: tsizes");
    for (j=0; j < ncols; j++) {
        tsizes[j] = csizes[cperm[j]]; // new column size after permuting
    }
    da_icopy(ncols, tsizes, csizes);
#else
    da_errexit("Invalid COLSORTBY criteria in defs.h.");
#endif

    /* set new column ids in docs */
    for (i=0; i < nnz; i++)
        rowind[i] = cperm[rowind[i]];

    da_csr_SortIndices(docs, DA_ROW);

    da_free((void**) &orderColNnz, &orderColVal, &cval, &tsizes, LTERM);

}
