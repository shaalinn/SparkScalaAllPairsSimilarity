/*!
 \file  pl2nn.c
 \brief This file contains parallel L2AP related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
void pl2nnReorderDocs(params_t * params, da_csr_t ** docs, idx_t * rsizes,
        idx_t * csizes, val_t * rwgts, val_t * cwgts, idx_t * rperm,
        idx_t * cperm);

/**
 * Main entry point to ParallelKnnIdxJoin.
 */
void pl2nnFindNeighbors(params_t * params)
{
    ssize_t i, j, k, nnz, nneighbs, mrlen, nmask, fsz, fnnz, innz, ntq;
    idx_t nrows, ncols, progressInd, pct, nqrows, ndrows, nqblks, ndblks, ninnz,
    qID, dID, tid;
    val_t v, b, bsq, simT;
    ptr_t *rowptr, *endptr = NULL, *blkptr = NULL;
    idx_t *rsizes, *csizes, *rperm, *cperm, *rowind;
    val_t *rowval, *rwgts = NULL, *rfiwgts = NULL, *cwgts = NULL, *ps = NULL,
            *rs = NULL, *prmax = NULL, *lens = NULL;
    da_idx_t **iidxs, *iidx, **fidxs, *fidx;
    da_tsearch_t **tsearch = NULL, *ts;
    da_csr_t *docs;

    /* pre-process the data */
    omp_set_num_threads(params->nthreads);
    preProcessData(params);
    da_csr_CompactColumns(params->docs);
    da_csr_PrintInfo(params->docs, "after compact: ", "\n");

    docs = params->docs;
    nrows = docs->nrows;		// num rows
    ncols = docs->ncols;		// num rows
    nnz = docs->rowptr[nrows];
    ninnz = da_min(params->ninnz, nnz);
    params->nqrows = da_min(nrows, params->nqrows);
    params->ndrows = da_min(nrows, params->ndrows);
    ndblks = ceil(nnz / (double) ninnz);	// max number of blocks
    simT = params->eps;
    ntq = 0;
    progressInd = 0;
#if L2APDP != L2APDP_HASH
    nmask = 0;
#endif

    /* allocate memory */
    da_startwctimer(params->timer_5);	// memory allocation time
    assert(nnz > 0);
    rsizes = da_imalloc(nrows, "pl2nnFindNeighbors: rsizes"); /** row sizes */
    csizes = da_ismalloc(ncols, 0.0, "pl2nnFindNeighbors: csizes"); /** col sizes */
    rperm = da_imalloc(nrows, "pl2nnFindNeighbors: rperm"); /** row permutation */
    cperm = da_imalloc(ncols, "pl2nnFindNeighbors: cperm"); /** col permutation */
    rwgts = da_vsmalloc(nrows, 0.0, "pl2nnFindNeighbors: rwgts"); /** max row value */
    rfiwgts = da_vsmalloc(nrows, 0.0, "pl2nnFindNeighbors: rfiwgts"); /** max row value in forward index */
    cwgts = da_vsmalloc(ncols, 0.0, "pl2nnFindNeighbors: cwgts"); /** max col value */
    ps = da_vmalloc(nrows, "pl2nnFindNeighbors: ps"); /** score prior to indexing threshold for each row */
    endptr = da_pmalloc(nrows, "pl2nnFindNeighbors: endptr"); /** where the inv index starts in each row */
    blkptr = da_pmalloc(ndblks + 1, "pl2nnFindNeighbors: blkptr"); /** pointer for the block starts */
#if PL2APRS == PL2APRS_STORE
    rs = da_vmalloc(nnz, "pl2nnFindNeighbors: rs"); /** remaining scores */
    lens = da_vmalloc(nnz, "pl2nnFindNeighbors: lens"); /** prefix lengths */
#else
    rs = da_vmalloc (nrows, "pl2nnFindNeighbors: rs"); /** first remaining score for each row */
#endif
#ifdef PL2AP_DP6
    prmax = da_vsmalloc(nnz, -1, "pl2nnFindNeighbors: prmax"); /** prefix max values for each row */
#endif
    tsearch = (da_tsearch_t **) da_malloc(
            sizeof(da_tsearch_t *) * params->nthreads,
            "ppl2nn2FindNeighbors: tsearch");
    da_stopwctimer(params->timer_5);	// memory allocation time

    // reorder matrix
    da_startwctimer(params->timer_3);// similarity search time -- includes indexing, cg and cv times
    pl2nnReorderDocs(params, &docs, rsizes, csizes, rwgts, cwgts, rperm, cperm);
    rowptr = docs->rowptr;

    /* Create the indexes and compute other meta-info in parallel */
    if (params->verbosity > 0)
    {
        printf("Indexing... ");
        fflush(stdout);
    }
    da_startwctimer(params->timer_4);	// indexing time
#pragma omp parallel
    {
        // find the index splits
#pragma omp for schedule(dynamic, NITERS)
        for (i = 0; i < nrows; ++i) {
            pl2nnFindIndexSplit(i, params, docs, endptr, rwgts, rfiwgts, cwgts,
                    ps, rs, lens, prmax);
        }
    }
    if (params->verbosity > 0)
    {
        printf("found splits... ");
        fflush(stdout);
    }

    /** find the number of blocks and store block pointers**/
    for (blkptr[0] = 0, i = 0, k = 0; k < nrows; i++)
    {
        for (j = 0; k < nrows;) {
            j += (rowptr[k + 1] - endptr[k]);
            if (j > ninnz)
                break;
            k++;
        }
        if (i >= ndblks)
            da_errexit("-ninnz parameter too small. Try a larger value.");
        blkptr[i + 1] = k;
    }
    ndblks = i;
    ASSERT(k == nrows);
    ASSERT(blkptr[i] == nrows);
    da_stopwctimer(params->timer_4);	// indexing time

    da_startwctimer(params->timer_5);	// memory allocation time
    iidxs = (da_idx_t **) da_malloc(ndblks * sizeof(da_idx_t *),
            "pl2nnFindNeighbors: iidxs");
    fidxs = (da_idx_t **) da_malloc(ndblks * sizeof(da_idx_t *),
            "pl2nnFindNeighbors: fidxs");
    for (i = 0; i < ndblks; ++i)
    {
        iidxs[i] = (da_idx_t *) da_malloc(sizeof(da_idx_t),
                "pl2nnFindNeighbors: iidxs[i]");
        memset(iidxs[i], 0, sizeof(da_idx_t));
        fidxs[i] = (da_idx_t *) da_malloc(sizeof(da_idx_t),
                "pl2nnFindNeighbors: iidxs[i]");
        memset(iidxs[i], 0, sizeof(da_idx_t));
    }
    da_stopwctimer(params->timer_5);	// memory allocation time

    da_startwctimer(params->timer_4);	// indexing time
#pragma omp parallel
    {
        for (i = 0; i < ndblks; i++) {
            pl2nnIndexRows(blkptr[i], blkptr[i + 1] - blkptr[i], params, docs,
                    endptr, iidxs[i], fidxs[i]);
        }
    }
    if (params->verbosity > 0)
    {
        printf("done.\n");
        fflush(stdout);
    }
    da_stopwctimer(params->timer_4);	// indexing time


#ifdef IDX_SIZE
    if (params->verbosity > 0)
    {
        printf ("done.\n");
        printf ("Index sizes: \n");
        fflush (stdout);
    }

    fsz = fnnz = innz = 0;
    for (i = 0; i < ndblks; i++)
    {
        iidx = iidxs[i];
        fidx = fidxs[i];
        printf ("Index %7zu\tFwdIdx: sz %zu nnz %zu\tInvIdx: sz %zu nnz %zu\n",
                i + 1, fidx->sz, fidx->ptr[fidx->sz], iidx->sz,
                iidx->ptr[iidx->sz]);
        fsz += fidx->sz;
        fnnz += fidx->ptr[fidx->sz];
        innz += iidx->ptr[iidx->sz];
        params->indexSize += iidx->ptr[iidx->sz];
    }
    v = i;
    printf ("Mean: nIdx %zu\tFwdIdx: sz %zu nnz %zu\tInvIdx: sz %zu nnz %zu\n",
            i, (size_t) (fsz / v), (size_t) (fnnz / v), iidxs[0]->sz,
            (size_t) (innz / v));
    da_stopwctimer (params->timer_3);// similarity search time -- includes indexing, cg and cv times
    /* free memory */
    // threads tmp data
    for (i = 0; i < params->nthreads; ++i)
    {
        ts = tsearch[i];
        da_free ((void **) &ts->accum, &ts->cands, &ts->hashlen, &ts->hashval,
                &ts->hashmax, &ts->hashtbl, &ts->hashdat, &ts->nbrslst, &ts,
                LTERM);
    }
    // indexes
    for (i = 0; i < ndblks; ++i)
    {
        da_free ((void **) &iidxs[i]->ptr, &iidxs[i]->ind, &iidxs[i]->len,
                &iidxs[i]->val, &iidxs[i], &fidxs[i]->ptr, &fidxs[i]->ind,
                &fidxs[i]->len, &fidxs[i]->val, &fidxs[i], LTERM);
    }
    da_free ((void **) &rsizes, &csizes, &rwgts, &rfiwgts, &cwgts, &endptr,
            &blkptr, &ps, &tsearch, &rs, &prmax, &iidxs, &fidxs, &lens, LTERM);
    return;
#endif



    // allocate working memory for each thread
    da_startwctimer(params->timer_5);	// memory allocation time
#if L2APDP != L2APDP_HASH
    mrlen = da_imax(nrows, rsizes);
#endif
    for (ndrows = 0, i = 0; i < ndblks; i++)
    {
        if (blkptr[i + 1] - blkptr[i] > ndrows)
            ndrows = blkptr[i + 1] - blkptr[i];
    }
    for (i = 0; i < params->nthreads; ++i)
    {
        tsearch[i] = ts = (da_tsearch_t *) da_malloc(sizeof(da_tsearch_t),
                "pl2nn2FindNeighbors: tsearch[i]");
        memset(ts, 0, sizeof(da_tsearch_t));
        ts->nrows = nrows;
        ts->ncols = ncols;
        ts->accum = da_asmalloc(ndrows, -1,
                "pl2nn2FindNeighbors: tsearch[i]->accum");
        ts->cands = da_ismalloc(ndrows, 0,
                "pl2nn2FindNeighbors: tsearch[i]->cands");
#ifdef PL2AP_SUP
        ts->start = da_pmalloc (ncols, "pl2nn2FindNeighbors: tsearch[i]->start");
#endif
        ts->hashval = da_vsmalloc(da_max(ncols, HTSIZE), 0,
                "pl2nn2FindNeighbors: tsearch[i]->hashval");
        ts->hashlen = da_vsmalloc(da_max(ncols, HTSIZE), 0,
                "pl2nn2FindNeighbors: tsearch[i]->hashlen");
        ts->hashmax = da_vsmalloc(da_max(ncols, HTSIZE), 0,
                "pl2nn2FindNeighbors: tsearch[i]->hashmax");
        ts->nbrsz = NINITNEIGHBORS * 1.25 * nrows / params->nthreads; /** initial neighbors space **/
        ts->nbrnnz = 0;
#ifndef NO_OUTPUT
        ts->nbrslst = da_smalloc(ts->nbrsz,
                "pl2nn2FindNeighbors: tsearch[i]->nbrslst");
#endif
#if L2APDP != L2APDP_HASH
        ts->hashdat = da_l2aphmalloc(HTSIZE + mrlen,
                "pl2nn2FindNeighbors: tsearch[i]->hashdat");
        ts->hashtbl = da_ismalloc(HTSIZE + mrlen, -1,
                "pl2nn2FindNeighbors: tsearch[i]->hashtbl");// hash table for masked hash
#endif
    }
    da_stopwctimer(params->timer_5);	// memory allocation time


    if (params->verbosity > 0)
    {
        printf("Progress Indicator: ");
        fflush(stdout);
        da_progress_init(pct, progressInd, (nrows / params->ndrows));
    }

    for (i = 0; i < ndblks; i++)
    {
        dID = blkptr[i];
        iidx = iidxs[i];
        fidx = fidxs[i];
#ifdef PL2AP_SUP
        for (j = 0; j < params->nthreads; ++j)
        {
            da_pcopy (iidx->sz, iidx->ptr, tsearch[j]->start); /** copies of index start for each thread **/
        }
#endif
        ntq += nrows - dID;

        for (qID = dID; qID < nrows; qID += params->nqrows) {
            nqrows = da_min(params->nqrows, nrows - qID);

#pragma omp parallel
            {
#pragma omp for private(tid, ts) schedule(dynamic, NITERS)
                for (j = 0; j < nqrows; ++j) {
                    tid = omp_get_thread_num();
                    ts = tsearch[tid];
#if L2APDP == L2APDP_MIX || L2APDP == L2APDP_MASK
                    printf(
                            "Note: Invalid L2APDP setting for pl2nn. Maked mode not implemented.\n");
                    pl2nnFindMatches(qID + j, qID, dID, simT, docs, rs, lens,
                            prmax, ps, rfiwgts, cwgts, rwgts[qID + j], rsizes,
                            iidx, fidx, ts, endptr);
#endif

                }
            }

        }
        if (params->verbosity > 0 && (dID / params->ndrows) % progressInd == 0)
            da_progress_advance(pct);

    }
    if (params->verbosity > 0)
    {
        da_progress_finalize(pct);
        printf("\n");
    }
    da_stopwctimer(params->timer_3);// similarity search time -- includes indexing, cg and cv times
#ifndef NO_OUTPUT
    params->neighbors = da_csr_CombineAndPermuteLT(nrows, params->nthreads,
            tsearch, params->rperm);
    if (params->verbosity > 0)
        for (i = 0; i < params->nthreads; ++i) {
            printf("neighbors_t%zu: ", i + 1);
            if (tsearch[i]->nbrs)
                da_csr_PrintInfo(tsearch[i]->nbrs, "", "\n");
            else
                printf("%zu\n", tsearch[i]->nbrnnz);
        }
    da_csr_PrintInfo(params->neighbors, "neighbors: ", "\n");
    if (params->fpout)
    {
        if (params->symout)
            da_addHigherNeighbors(params->neighbors);
        if (params->verbosity > 0)
            printf("Writing neighborhood matrix to %s.\n", params->oFile);
        da_csr_Write(params->neighbors, params->oFile, params->fmtWrite, 1, 1);
    }
#endif


    /** finalize sim search **/
    // transfer counts
    for (i = 0; i < params->nthreads; ++i)
    {
        ts = tsearch[i];
        params->nSimPairs += ts->nSimPairs;
        params->nDotProducts += ts->nDotProducts;
        params->nDotProducts1 += ts->nDotProducts1;
        params->nCandidates += ts->nCandidates;
        nmask += ts->nmask;
#ifdef EXTRACOUNTS
        params->nPruneMinsize += ts->nPruneMinsize;
        params->nPruneDotP += ts->nPruneDotP;
        params->nPruneDotP2 += ts->nPruneDotP2;
        params->nPrunePscore += ts->nPrunePscore;
        params->nPruneLength += ts->nPruneLength;
        params->nPruneLength2 += ts->nPruneLength2;
#ifdef HASHCOUNTS
        params->nCacheHit += ts->nCacheHit;
        params->nCacheMiss += ts->nCacheMiss;
#endif
#ifdef ACCUMCOUNTS
        params->percAccum += ts->percAccum;
#endif
#endif
#ifdef EXTRATIMES
        if (params->verbosity > 0)
            printf("Thread %zu\tcg %.2f\tcv %.2f\n", i + 1, ts->timer_1,
                    ts->timer_2);
        if (ts->timer_1 > params->timer_1)
            params->timer_1 = ts->timer_1;
        if (ts->timer_2 > params->timer_2)
            params->timer_2 = ts->timer_2;
#endif

    }

#ifdef EXTRACOUNTS
    fsz = fnnz = innz = 0;
    for (i = 0; i < ndblks; i++)
    {
        iidx = iidxs[i];
        fidx = fidxs[i];
        printf ("Index %7zu\tFwdIdx: sz %zu nnz %zu\tInvIdx: sz %zu nnz %zu\n",
                i + 1, fidx->sz, fidx->ptr[fidx->sz], iidx->sz,
                iidx->ptr[iidx->sz]);
        fsz += fidx->sz;
        fnnz += fidx->ptr[fidx->sz];
        innz += iidx->ptr[iidx->sz];
        params->indexSize += iidx->ptr[iidx->sz];
    }
    v = i;
    printf ("Mean: nIdx %zu\tFwdIdx: sz %zu nnz %zu\tInvIdx: sz %zu nnz %zu\n",
            i, (size_t) (fsz / v), (size_t) (fnnz / v), iidxs[0]->sz,
            (size_t) (innz / v));
#ifdef ACCUMCOUNTS
    if (params->verbosity > 0)
    {
        printf ("Mean percent computed accumulation (no dp): %.6f\n",
                params->percAccum / params->nCandidates);
    }
    params->percAccum += (params->nDotProducts + params->nDotProducts1);// had 100% accumulation
    if (params->verbosity > 0)
    {
        printf ("Mean percent computed accumulation: %.6f\n",
                params->percAccum / params->nCandidates);
    }
#endif
#endif

    if (params->verbosity > 0)
    {
        printf("Number of executed queries: %zu\n", ntq);
    }

#if L2APDP != L2APDP_HASH
    if (params->verbosity > 0)
    {
        printf("Number of masked hash queries: %zu\n", nmask);
        printf("Percent rows using dense hash: %.2f\n",
                (ntq - nmask) / (double) (ntq));
        printf("Percent rows using masked hash: %.2f\n",
                nmask / (double) (ntq));
#ifdef HASHCOUNTS
        v = (params->nCacheHit + params->nCacheMiss);
        if (v == 0)
        v = 1;
        printf ("Masked hashtable hits:   %zu (count), %.2f (pcnt)\n",
                params->nCacheHit, params->nCacheHit / v);
        printf ("Masked hashtable misses: %zu (count), %.2f (pcnt)\n",
                params->nCacheMiss, params->nCacheMiss / v);
#endif
    }
#endif

    /* free memory */
    // threads tmp data
    for (i = 0; i < params->nthreads; ++i)
    {
        ts = tsearch[i];
        da_free((void **) &ts->accum, &ts->cands, &ts->hashlen, &ts->hashval,
                &ts->start, &ts->hashmax, &ts->hashtbl, &ts->hashdat,
                &ts->nbrslst, &ts, LTERM);
    }
    // indexes
    for (i = 0; i < ndblks; ++i)
    {
        da_free((void **) &iidxs[i]->ptr, &iidxs[i]->ind, &iidxs[i]->len,
                &iidxs[i]->val, &iidxs[i], &fidxs[i]->ptr, &fidxs[i]->ind,
                &fidxs[i]->len, &fidxs[i]->val, &fidxs[i], LTERM);
    }
    da_free((void **) &rsizes, &csizes, &rwgts, &rfiwgts, &cwgts, &endptr,
            &blkptr, &ps, &tsearch, &rs, &prmax, &iidxs, &fidxs, &lens, LTERM);
}




/**
 * Find matches for row id rid
 * rid = id for row in pmat, i.e. block of query rows
 */
void pl2nnFindMatches(idx_t rid, idx_t qstart, idx_t dstart, val_t simT,
        da_csr_t * docs, val_t * rs, val_t * lens, val_t * prmax, val_t * ps,
        val_t * rfiwgts, val_t * cwgts, val_t rwgt, idx_t * rsizes,
        da_idx_t * iidx, da_idx_t * fidx, da_tsearch_t * ts, ptr_t * endptr)
{

    da_errexit("This mode has not yet been implemented.");

#ifdef EXTRATIMES
    da_startwctimer(ts->timer_1);	// find candidates time
#endif
    ssize_t i, j, b;
    idx_t k, ncands, cid, lcid;
    ptr_t sptr, eptr;
    ptr_t *matptr, *idxptr, *starts;
    idx_t *matind, *idxind, *cands;
    val_t simTe, myval, cval, add, bnd, rs1, rs2, bsq, acc;
    val_t *matval, *idxval, *hashval, *hashlen, *hashmax, *idxlen;
    accum_t *accum;

    matptr = docs->rowptr;	// index in rowind/rowval where each row starts
    matind = docs->rowind;	// col id for each nnz
    matval = docs->rowval;	// val for each nnz
    idxind = iidx->ind;
    idxval = iidx->val;
    idxlen = iidx->len;
    idxptr = iidx->ptr;		// starting points for lists in the inv index
#ifdef PL2AP_MINSZ
#ifdef PL2AP_SUP
    starts = ts->start;		// starting points for lists in the inv index
#else
    starts = iidx->ptr;		// starting points for lists in the inv index
#endif
#endif
    accum = ts->accum;		// accumulator array
    hashval = ts->hashval;	// dense version of query vector
    hashlen = ts->hashlen;	// dense version of query vector suffix lengths
    hashmax = ts->hashmax;	// dense version of query vector prefix max values
    cands = ts->cands;		// candidates array
    sptr = matptr[rid];
    eptr = matptr[rid + 1];

    ncands = 0;			// number of candidates for this doc
    simTe = simT - LENBNDEPS;
#if PL2APRS == PL2APRS_LIVE
    rs1 = rs[rid];
#else
    rs1 = eptr == sptr ? 0.0 : rs[eptr - 1];
#endif
#if PL2APRS != PL2APRS_STORE
    bsq = rs2 = 1.0;
#endif

#ifdef PL2AP_MINSZ
    const accum_t minsize = ceil((simT / rwgt) * (simT / rwgt));
    if (rs1 < simT || minsize > docs->ncols)
    {
#else
        if (rs1 < eps)
        {
#endif
#ifdef EXTRATIMES
        da_stopwctimer(ts->timer_1);	// find candidates time
#endif
#ifdef EXTRACOUNTS
        ts->nPruneMinsize++;
#endif
        return;
    }

    for (i = eptr - 1; i >= sptr; i--) {
        k = matind[i];
        myval = matval[i];
        hashval[k] = myval;

#if PL2APRS == PL2APRS_STORE
        bnd = rs[i];
        bsq = hashlen[k] = lens[i];
#else
        bnd = da_min (rs1, bsq);
        rs2 -= (double) myval *myval;
        bsq = sqrt (rs2);
        hashlen[k] = bsq;
#endif
#ifdef PL2AP_DP6
        hashmax[k] = prmax[i];
#endif

#ifdef PL2AP_MINSZ
        const ptr_t b = starts[k];
        for (j = b; j < idxptr[k + 1]; j++) {
            cid = idxind[j];
            if (cid >= rid || rsizes[cid] >= minsize) {
#ifdef EXTRACOUNTS
                /**
                 * Note that, since query vector is parsed multiple times in pl2nn/pl2nn2,
                 * nPruneMinsize will differ (be larger than that of pl2nn1), unless only one
                 * index is created.
                 */
                ts->nPruneMinsize++;
#endif
                break;
            }
        }
#ifdef PL2AP_SUP
        if (j > b)
        {
            starts[k] += j - b;
        }
#endif
#else
        j = idxptr[k];
#endif

        for (; j < idxptr[k + 1]; j++) {
            cid = idxind[j];	// potential candidate from the inv index
            if (cid >= rid) {
                break;
            }
            lcid = cid - dstart;
            acc = accum[lcid];

            if (acc > 0.0) {
                // accumulate
                accum[lcid] = acc += myval * idxval[j];
                // check l2-norm bound
                if (acc + bsq * idxlen[j] < simTe) {
                    accum[lcid] = -2;
#ifdef EXTRACOUNTS
                    ts->nPruneLength++;
#ifdef ACCUMCOUNTS
                    ts->percAccum +=
                    computePercAccum (rid, cid, k, -1, docs, endptr);
#endif
#endif
                }
            } else if (bnd >= simT && acc > -2) {
                // accumulate
                if (acc == -1) {
                    accum[lcid] = acc = myval * idxval[j];
                    cands[ncands++] = cid - dstart;
                } else
                    accum[lcid] = acc += myval * idxval[j];
                // check l2-norm bound
                if (acc + bsq * idxlen[j] < simTe) {
                    accum[lcid] = -2;
#ifdef EXTRACOUNTS
                    ts->nPruneLength++;
#ifdef ACCUMCOUNTS
                    ts->percAccum +=
                    computePercAccum (rid, cid, k, -1, docs, endptr);
#endif
#endif
                }
            }

        }

#if PL2APRS == PL2APRS_LIVE
        rs1 -= myval * cwgts[k];
#endif

    }
#ifdef EXTRATIMES
    da_stopwctimer(ts->timer_1);	// find candidates time

    da_startwctimer(ts->timer_2);	// process candidates time
#endif
    ts->nCandidates += ncands;
    ts->nDotProducts += pl2nnProcessCandidates(rid, qstart, dstart, simT,
            eptr - sptr, rwgt, rfiwgts, ps, hashval, hashlen, hashmax, fidx,
            ncands, cands, accum, ts, docs, endptr);
#ifdef EXTRATIMES
    da_stopwctimer(ts->timer_2);	// process candidates time

    da_startwctimer(ts->timer_1);	// find candidates time
#endif
    // reset hashval
    for (i = sptr; i < eptr; i++)
        hashval[matind[i]] = 0;

#ifdef EXTRATIMES
    da_stopwctimer(ts->timer_1);	// find candidates time
#endif
}




/**
 * Process identified candidates
 */
idx_t
pl2nnProcessCandidates(idx_t rid, idx_t qstart, idx_t dstart, val_t simT,
        size_t sz, val_t rwgt, val_t * rfiwgts, val_t * ps, val_t * hashval,
        val_t * hashlen, val_t * hashmax, da_idx_t * fidx, idx_t ncands,
        idx_t * cands, accum_t * accum, da_tsearch_t * ts, da_csr_t * docs,
        ptr_t * endptr)
{

    ssize_t i, j, k, p, r, ls, cid, lcid, dps = 0, nbrnnz;
    val_t crfiwgt;
    double cval, simTe;
    ptr_t *dptr;
    idx_t *dind;
    val_t *dval, *dlen;

    simTe = simT - LENBNDEPS;
    dptr = fidx->ptr;
    dind = fidx->ind;
    dval = fidx->val;
    dlen = fidx->len;

    // grow neighbors list if necessary
#ifndef NO_OUTPUT
    if (ts->nbrnnz + ncands > ts->nbrsz)
    {
        ts->nbrsz += da_max(ncands, ts->nbrnnz / 4);
        ts->nbrslst = da_srealloc(ts->nbrslst, ts->nbrsz,
                "pl2nnProcessCandidates: ts->nbrslst");
    }
#endif

    for (i = 0; i < ncands; i++)
    {
        lcid = cands[i];

        if (accum[lcid] <= 0.0) {
            accum[lcid] = -1;
            continue;
        }

        cval = accum[lcid];
        accum[lcid] = -1;
        cid = lcid + dstart;

        if (cval + ((double) ps[cid]) < simT) {
#ifdef EXTRACOUNTS
            ts->nPrunePscore++;
#ifdef ACCUMCOUNTS
            ts->percAccum += computePercAccum (rid, cid, 0, -1, docs, endptr);
#endif
#endif
            continue;
        }

        crfiwgt = rfiwgts[cid];
        if (cval
                + ((double) da_min(dptr[lcid + 1] - dptr[lcid], sz) * crfiwgt
                        * rwgt) < simT) {
#ifdef EXTRACOUNTS
            ts->nPruneDotP++;
#ifdef ACCUMCOUNTS
            ts->percAccum += computePercAccum (rid, cid, 0, -1, docs, endptr);
#endif
#endif
            continue;
        }

        r = dptr[lcid];

        // advance to next term in common
        for (j = dptr[lcid + 1] - 1; j >= r; j--)
            if (hashval[dind[j]] > 0.0)
                break;

        if (j < r) {
            ts->nDotProducts1++; /** full dot product before verification **/
            goto check;
        }

#ifdef PL2AP_DP6
        if (cval + ((double) da_min(j - r + 1, sz) * crfiwgt * hashmax[dind[j]])
                < simT) {
#ifdef EXTRACOUNTS
            ts->nPruneDotP2++;
#ifdef ACCUMCOUNTS
            ts->percAccum +=
            computePercAccum (rid, cid, 0, dind[j], docs, endptr);
#endif
#endif
            continue;
        }
#endif

        for (; j >= r; j--)
            if (hashval[dind[j]] > 0.0) {
                cval += (double) hashval[dind[j]] * dval[j];
                if (j > r
                        && cval + ((double) hashlen[dind[j]] * dlen[j])
                                < simTe) {
#ifdef EXTRACOUNTS
                    ts->nPruneLength2++;
#ifdef ACCUMCOUNTS
                    ts->percAccum +=
                    computePercAccum (rid, cid, 0, dind[j], docs, endptr);
#endif
#endif
                    goto nexcid;
                }
            }

        dps++;

        check:			// check similarity value

        if (cval >= simT) {
            ts->nSimPairs++;
            //add to local neighbors list
#ifndef NO_OUTPUT
            ts->nbrslst[ts->nbrnnz].i = rid;
            ts->nbrslst[ts->nbrnnz].j = cid;
            ts->nbrslst[ts->nbrnnz++].val = cval;
#endif
        }
        nexcid: continue;

    }

    return dps;


}





/**
 * Find where the forward index stops in the row and compute other necessary scores
 */
void pl2nnFindIndexSplit(idx_t rid, params_t * params, da_csr_t * docs,
        ptr_t * endptr, val_t * rwgts, val_t * rfiwgts, val_t * cwgts,
        val_t * ps, val_t * rs, val_t * lens, val_t * prmax)
{
    ssize_t j, k;
    val_t myval, maxrval, maxcval;
    double b1, b2, b2l, sqb2, rs1, pscore, simT;
    ptr_t *rowptr;
    idx_t *rowind;
    val_t *rowval;

    rowptr = docs->rowptr;	// index in rowind/rowval where each row starts
    rowind = docs->rowind;	// colid for each nnz
    rowval = docs->rowval;	// val for each nnz
    simT = params->eps;
    pscore = 0.0;
    myval = b1 = maxrval = b2 = sqb2 = b2l = rs1 = 0.0;

    // index terms in row i
    for (endptr[rid] = rowptr[rid + 1], j = rowptr[rid]; j < rowptr[rid + 1];
            j++) {
#if PL2APRS == PL2APRS_STORE
        lens[j] = sqb2;
#endif
        myval = rowval[j];
        maxcval = cwgts[rowind[j]];
        b1 += (double) myval * da_min(rwgts[rid], maxcval);
        // adjust b2 estimate of prefix similarity with all other vectors
        b2 += (double) myval * myval;
        sqb2 = sqrt(b2);

        // check whether to start indexing
        if (da_min (b1, sqb2) >= simT)
            break;

        // compute remaining score
        rs1 += myval * maxcval;
#if PL2APRS == PL2APRS_STORE
        rs[j] = da_min(rs1, sqb2);
#endif
        //update pscore
        pscore = da_min(b1, sqb2);
        b2l = b2;
        if (myval > maxrval)
            maxrval = myval;
#ifdef PL2AP_DP6
        prmax[j] = maxrval;
#endif
    }
    // truncate the rest of the vector, since we're going to index it.
    endptr[rid] = j;
    // update max forward index row val for row rid
    rfiwgts[rid] = maxrval;
    //store the pscore for this row
    ps[rid] = pscore;

    b2 -= (double) myval * myval;
    sqb2 = b2 > 0 ? sqrt(b2) : 0.0;
    // continue storing prefix max values and computing rs1
    for (; j < rowptr[rid + 1]; j++)
    {
#if PL2APRS == PL2APRS_STORE
        lens[j] = sqb2;
#endif
        myval = rowval[j];
        rs1 += myval * cwgts[rowind[j]];
#if PL2APRS == PL2APRS_STORE
        b2 += (double) myval * myval;
        sqb2 = sqrt(b2);
        rs[j] = da_min(rs1, sqb2);
#endif
#ifdef PL2AP_DP6
        if (myval > maxrval)
            maxrval = myval;
        prmax[j] = maxrval;
#endif
    }
#if PL2APRS == PL2APRS_LIVE
    // store initial remaining score for each row
    rs[rid] = rs1;
#endif
}



/**
 * Create forward and inverse indexes for a block of rows
 */
void pl2nnIndexRows(idx_t start, idx_t nrows, params_t * params,
        da_csr_t * docs,
        ptr_t * endptr, da_idx_t * iidx, da_idx_t * fidx)
{
    idx_t rid, end;
    ssize_t i, j, k, n, innz, fnnz, ncols;
    double b2, sqb2;
    ptr_t *rowptr, *iptr, *fptr;
    idx_t *rowind, *iind, *find;
    val_t *rowval, *ival, *ilen, *fval, *flen;

    end = start + nrows;
    ncols = docs->ncols;
    rowptr = docs->rowptr;	// index in rowind/rowval where each row starts
    rowind = docs->rowind;	// colid for each nnz
    rowval = docs->rowval;	// val for each nnz

    iidx->ptr = iptr = da_psmalloc(ncols + 1, 0,
            "pl2nnFindNeighbors: invIdx->ptr");
    for (innz = 0, fnnz = 0, rid = start; rid < end; rid++)
    {
        fnnz += endptr[rid] - rowptr[rid];
        innz += rowptr[rid + 1] - endptr[rid];
        for (j = endptr[rid]; j < rowptr[rid + 1]; j++) {
            iptr[rowind[j]]++;
        }
    }

    CSRMAKE(i, ncols, iptr);
    iidx->sz = ncols;
    iidx->ind = iind = da_ismalloc(innz, 0, "pl2nnFindNeighbors: iidx->ind");
    iidx->val = ival = da_vsmalloc(innz, 0.0, "pl2nnFindNeighbors: iidx->val");
    iidx->len = ilen = da_vsmalloc(innz, 0.0, "pl2nnFindNeighbors: iidx->len");
    fidx->sz = nrows;
    fidx->ptr = fptr = da_psmalloc(nrows + 1, 0,
            "pl2nnFindNeighbors: fidx->ptr");
    fidx->ind = find = da_ismalloc(fnnz, 0, "pl2nnFindNeighbors: fidx->ind");
    fidx->val = fval = da_vsmalloc(fnnz, 0.0, "pl2nnFindNeighbors: fidx->val");
    fidx->len = flen = da_vsmalloc(fnnz, 0.0, "pl2nnFindNeighbors: fidx->len");

    for (fnnz = 0, rid = start; rid < end; rid++)
    {
        n = endptr[rid] - rowptr[rid];
        da_icopy(n, rowind + rowptr[rid], fidx->ind + fnnz);
        da_vcopy(n, rowval + rowptr[rid], fidx->val + fnnz);
        for (b2 = sqb2 = 0.0, j = rowptr[rid]; j < endptr[rid]; j++, fnnz++) {
            flen[fnnz] = sqb2;
            b2 += (double) rowval[j] * rowval[j];
            sqb2 = sqrt(b2);
        }
        fptr[rid - start + 1] = fnnz;

        // index the remaining part of the vector
        for (j = endptr[rid]; j < rowptr[rid + 1]; j++) {
            k = rowind[j];
            ilen[iptr[k]] = sqb2;
            b2 += (double) rowval[j] * rowval[j];
            sqb2 = sqrt(b2);
            iind[iptr[k]] = rid;
            ival[iptr[k]] = rowval[j];
            iptr[k]++;
        }

    }

    CSRSHIFT(i, ncols, iptr);

}

/**
 * Reorder the document matrix in decreasing number of column nnz
 * and decreasing order of row max val. Also compute and store
 * necessary stats about the data
 * Store initial order in rperm/cperm.
 */
void pl2nnReorderDocs(params_t * params, da_csr_t ** docs, idx_t * rsizes,
        idx_t * csizes, val_t * rwgts, val_t * cwgts, idx_t * rperm,
        idx_t * cperm) {
    ssize_t i, j, k, nnz;
    idx_t nrows, ncols;
    da_iikv_t *orderColNnz = NULL;
    da_ivkv_t *orderMaxVal = NULL;
    char pv = 0;
    ptr_t *rowptr, *nrowptr;
    idx_t *rowind, *nrowind;
    val_t v, *rowval, *nrowval;

    nrows = (*docs)->nrows;
    ncols = (*docs)->ncols;
    rowptr = (*docs)->rowptr;
    rowind = (*docs)->rowind;
    rowval = (*docs)->rowval;
    nnz = rowptr[nrows];

    //record col nnzs
    da_iset(ncols, 0, csizes);
    for (i = 0; i < nnz; i++)
        csizes[rowind[i]]++;

    //assign memory for ordering
    orderMaxVal = da_ivkvmalloc(nrows, "l2apReorderDocs: orderMaxVal");
    nrowptr = da_pmalloc(nrows + 1, "l2apReorderDocs: nrowptr");
    nrowind = da_imalloc(nnz, "l2apReorderDocs: nrowind");
    nrowval = da_vmalloc(nnz, "l2apReorderDocs: nrowval");

    orderColNnz = da_iikvmalloc(ncols, "l2apReorderDocs: orderColNnz");
    //get new column order
    for (j = 0; j < ncols; j++) {
        orderColNnz[j].key = j;
        orderColNnz[j].val = csizes[j];
    }
    da_iikvsortd(ncols, orderColNnz);   // sort columns
    for (j = 0; j < ncols; j++) {
        cperm[orderColNnz[j].key] = j;  // store column permutation
        csizes[j] = orderColNnz[j].val; // new column size after permuting
    }

    // set new column ids in docs
    for (i = 0; i < nnz; i++)
        rowind[i] = cperm[rowind[i]];
    // get max weights for rows and columns + row sizes
    for (i = 0; i < nrows; i++) {
        rsizes[i] = (idx_t) (rowptr[i + 1] - rowptr[i]);    //record row nnzs
        for (j = rowptr[i]; j < rowptr[i + 1]; j++) {
            v = rowval[j];
            if (v > rwgts[i])
                rwgts[i] = v;   // store max val for row
            if (v > cwgts[rowind[j]])
                cwgts[rowind[j]] = v;   //store max val for col
        }
        //set up kv array to get new row order
        orderMaxVal[i].key = i;
        orderMaxVal[i].val = rwgts[i];
    }
    da_ivkvsortd(nrows, orderMaxVal);   // sort rows
    for (i = 0; i < nrows; i++) {
        rperm[i] = orderMaxVal[i].key;  // store new order
        rwgts[i] = orderMaxVal[i].val;  // new row weight after permuting
    }
    // permute rows in matrix
    for (nrowptr[0] = 0, nnz = 0, j = 0, i = 0; i < nrows; i++) {
        k = rperm[i];
        da_icopy(rowptr[k + 1] - rowptr[k], rowind + rowptr[k], nrowind + nnz);
        da_vcopy(rowptr[k + 1] - rowptr[k], rowval + rowptr[k], nrowval + nnz);
        rsizes[i] = rowptr[k + 1] - rowptr[k];
        nnz += rsizes[i];
        nrowptr[++j] = nnz;
    }
    da_free((void **) &(*docs)->rowptr, &(*docs)->rowind, &(*docs)->rowval,
            LTERM);
    (*docs)->rowptr = nrowptr;
    (*docs)->rowind = nrowind;
    (*docs)->rowval = nrowval;

    da_csr_SortIndices((*docs), DA_ROW);

    da_free((void **) &orderColNnz, &orderMaxVal, &cperm, LTERM);

    params->rperm = rperm;

}
