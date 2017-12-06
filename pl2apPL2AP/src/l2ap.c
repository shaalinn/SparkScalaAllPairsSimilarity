/*!
 \file  l2ap.c
 \brief This file contains l2ap related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
idx_t l2apProcessCandidates (params_t * params, idx_t rid, da_csr_t * docs,
        val_t * lengths, val_t * sums, ptr_t * endptr, idx_t ncands,
        idx_t * cands, accum_t * accum, val_t * hashval, val_t * hashwgt,
        val_t * hashlen, val_t * hashsum, idx_t * hashsz, val_t * rwgts,
        val_t * rfiwgts, val_t * ps);
void l2apIndexRow (idx_t rid, params_t * params, da_csr_t * docs,
        ptr_t * endptr,
        da_invIdxJ_t * invIdx, val_t * rwgts, val_t * rfiwgts, val_t * cwgts,
        val_t * criwgts, val_t * hashwgt, val_t * ps);
void l2apFindMatches (idx_t rid, params_t * params, idx_t * cands,
        da_invIdxJ_t * idx, val_t * lengths, val_t * sums, da_csr_t * docs,
        idx_t * rsizes, val_t * rwgts, val_t * rfiwgts, val_t * cwgts,
        val_t * criwgts, val_t * ps, ptr_t * endptr, val_t * hashval,
        val_t * hashwgt, val_t * hashlen, val_t * hashsum, idx_t * hashsz);
void l2apReorderDocs (params_t * params, da_csr_t ** docs, idx_t * rsizes,
        idx_t * csizes, val_t * rwgts, val_t * cwgts, idx_t * rperm,
        idx_t * cperm);

/**
 * Main entry point to L2AP.
 */
void
l2apFindNeighbors (params_t * params)
{

    ssize_t i, j, k, nnz;
    idx_t nrows, ncols, mrlen, mclen, progressInd;
    ptr_t *rowptr, *endptr;
    idx_t *cands, *rsizes, *csizes, *rperm, *cperm, *hashsz = NULL;
    val_t maxrval, myval;
    val_t *rwgts = NULL, *rfiwgts = NULL, *hashwgt = NULL, *cwgts = NULL,
            *criwgts = NULL, *hashval = NULL, *hashlen = NULL, *hashsum = NULL,
            *ps = NULL, *lengths = NULL, *sums = NULL;
    da_invIdxJ_t *invIdx = NULL;//alternative inverted index for value-driven similarities//
    da_csr_t *docs = params->docs;

    // allocate memory
    params->nim = 1; /** default to storing neighbors in memory */
    simSearchSetup(params);
    nrows = docs->nrows;		// num rows
    ncols = docs->ncols;		// num cols
    params->eps2 = params->eps * params->eps;
    progressInd = ceil(nrows / (float) 20);

    rsizes = da_imalloc(nrows, "l2apFindNeighbors: rsizes");
    csizes = da_ismalloc(ncols, 0.0, "l2apFindNeighbors: csizes");
    rperm = da_imalloc(nrows, "l2apFindNeighbors: rperm");
    cperm = da_imalloc(ncols, "l2apFindNeighbors: cperm");
    rwgts = da_vsmalloc(nrows, 0.0, "l2apFindNeighbors: rwgts");
    rfiwgts = da_vmalloc(nrows, "l2apFindNeighbors: rfiwgts");
    cwgts = da_vsmalloc(ncols, 0.0, "l2apFindNeighbors: cwgts");
#if defined(RS1) && defined(RS3)
    criwgts = da_vsmalloc (ncols, 0.0, "l2apFindNeighbors: criwgts");
#endif
    hashval = da_vsmalloc(ncols, 0, "l2apFindNeighbors: hashval");// hash values of x
#if defined(DP7) || defined(DP8)
            hashsz = da_imalloc (ncols, "l2apFindNeighbors: hashsz");
#endif
#if defined(L2CV) || defined(LENCV)
    hashlen = da_vsmalloc(ncols, 0, "l2apFindNeighbors: hashlen");// hash suffix lengths of x
#endif
#if defined(DP3) || defined(DP4)
    hashsum = da_vsmalloc (ncols, 0, "l2apFindNeighbors: hashsum");	// hash suffix sums of x
#endif
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
    hashwgt = da_vmalloc(ncols, "l2apFindNeighbors: hashwgt");
#endif
#ifdef PSCV
    ps = da_vmalloc(nrows, "l2ap2FindNeighbors: ps");// score prior to indexing threshold for each row
#endif

    // reorder matrix
    l2apReorderDocs(params, &docs, rsizes, csizes, rwgts, cwgts, rperm, cperm);

    nnz = docs->rowptr[nrows];	// nnz in docs
    rowptr = docs->rowptr;	// index in rowind/rowval where each row starts

    // allocate memory for candidates and the inverted indexes -- needs nnz & col counts
    cands = da_imalloc(nrows, "l2apFindNeighbors: candidates");
#if defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4)
    sums = da_vsmalloc (nnz, 0.0, "l2apFindNeighbors: sums");// keep track of suffix sums
#endif
    lengths = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: lengths");// keep track of suffix lengths
    endptr = da_pmalloc(nrows, "l2apFindNeighbors: endptr");
    invIdx = da_malloc(sizeof(da_invIdxJ_t), "l2apFindNeighbors: invIdx");
    invIdx->ids = da_ismalloc(nnz, 0, "l2apFindNeighbors: invIdx->ids");
    invIdx->vals = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: invIdx->ids");
    invIdx->lens = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: invIdx->ids");
    invIdx->starts = da_pmalloc(ncols, "l2apFindNeighbors: invIdx->starts");
    invIdx->ends = da_pmalloc(ncols, "l2apFindNeighbors: invIdx->ends");
    for (j = 0, k = 0; j < ncols; j++)
    {
        invIdx->starts[j] = invIdx->ends[j] = k;
        k += csizes[j];
    }
    invIdx->accum = da_asmalloc(nrows, -1, "l2apFindNeighbors: invIdx->accum");

#ifdef EXTRACOUNTS
    // find the average length of rows & columns
    mrlen = ceil (da_isum (nrows, rsizes, 1) / (val_t) nrows);
    mclen = ceil (da_isum (ncols, csizes, 1) / (val_t) ncols);
    if (params->verbosity > 0)
    printf ("mrlen: %d, mclen: %d.\n", mrlen, mclen);
#endif

    if (params->verbosity > 0)
        printCompileChoices();

    // for each x \in V do
    da_startwctimer(params->timer_3);// similarity search time -- includes indexing, cg and cv times
    if (params->verbosity > 0)
        printf("Progress Indicator (no. of records): ");
    for (i = 0; i < nrows; i++)
    {
        endptr[i] = rowptr[i + 1];
        // find matches for row i
        l2apFindMatches(i, params, cands, invIdx, lengths, sums, docs, rsizes,
                rwgts, rfiwgts, cwgts, criwgts, ps, endptr, hashval, hashwgt,
                hashlen, hashsum, hashsz);

#ifdef EXTRATIMES
        da_startwctimer(params->timer_4);	// indexing time
        l2apIndexRow(i, params, docs, endptr, invIdx, rwgts, rfiwgts, cwgts,
                criwgts, hashwgt, ps);
        da_stopwctimer(params->timer_4);	// indexing time
#else
                l2apIndexRow (i, params, docs, endptr, invIdx, rwgts, rfiwgts, cwgts,
                        criwgts, hashwgt, ps);
#endif

        if (params->verbosity > 0 && i % progressInd == 0) {
            printf("%d%%..", (int) (100 * i / (float) nrows));
            fflush(stdout);
        }
    }
    if (params->verbosity > 0)
        printf("100%%\n");
#ifdef EXTRACOUNTS
    for (j = 0; j < ncols; j++)
    params->indexSize += invIdx->ends[j] - invIdx->starts[j];
#endif
    da_stopwctimer(params->timer_3);// similarity search time -- includes indexing, cg and cv times

    simSearchFinalize(params);
    da_free((void **) &rsizes, &csizes, &rwgts, &rfiwgts, &cwgts, &criwgts,
            &cands, &endptr, &hashval, &hashwgt, &hashlen, &hashsum, &hashsz,
            &lengths, &sums, &ps, &invIdx->ids, &invIdx->vals, &invIdx->ends,
            &invIdx->starts, &invIdx->accum, &invIdx->lens, &invIdx, LTERM);

}



/**
 * Find matches for row id rid
 */
void
l2apFindMatches (idx_t rid, params_t * params, idx_t * cands,
        da_invIdxJ_t * idx, val_t * lengths, val_t * sums, da_csr_t * docs,
        idx_t * rsizes, val_t * rwgts, val_t * rfiwgts, val_t * cwgts,
        val_t * criwgts, val_t * ps, ptr_t * endptr, val_t * hashval,
        val_t * hashwgt, val_t * hashlen, val_t * hashsum, idx_t * hashsz)
{

#ifdef EXTRATIMES
    da_startwctimer(params->timer_1);	// find candidates time
#endif
    ssize_t i, j, s, ls;
    idx_t k, ncands, cid;
    ptr_t *rowptr, *starts, *ends;
    idx_t *rowind, *idxids;
    val_t v, minsize, myval;
    val_t *rowval, *idxvals, *idxlens;
    accum_t *accum;
    double simT, simTe, cval, sum, bnd, rs1, rs2, bsq, acc;

    rowptr = docs->rowptr;	// index in rowind/rowval where each row starts
    rowind = docs->rowind;	// col id for each nnz
    rowval = docs->rowval;	// val for each nnz
    idxids = idx->ids;
    idxvals = idx->vals;
    starts = idx->starts;		// starting points for lists in the inv index
    ends = idx->ends;		// ending points for lists in the inv index
    accum = idx->accum;		// accumulator array
    ncands = 0;			// number of candidates for this doc
    simT = params->eps;
#if defined(IDXL2)
    rs2 = 1.0;
    bsq = 1.0;
#elif defined(IDXLEN)
    rs2 = 1.0;
#endif
#if defined(L2CG) || defined(LENCG)
    idxlens = idx->lens;
#endif
#if defined(L2CG)
    simTe = params->eps - LENBNDEPS;
#endif
#ifdef RS1
    rs1 = 0.0;
#endif

    // define vars for hashing DP data
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
    v = 0.0;
#endif
#if defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4)
    sum = 0.0;
#endif
#if defined(DP7) || defined(DP8)
    k = 1;
#endif


    // compute rs1, the remaining score
#if defined(RS1) || defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4) \
        || defined(DP6) || defined(DP7) || defined(DP8)
    for (i = rowptr[rid]; i < rowptr[rid + 1]; i++)
    {
#ifdef RS1
#ifdef RS3
        rs1 += rowval[i] * criwgts[rowind[i]];	// val * max col weight
#else
        rs1 += rowval[i] * cwgts[rowind[i]];	// val * max col weight
#endif
#endif
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
        if (rowval[i] > v)
            v = rowval[i];
        hashwgt[rowind[i]] = v;	// prefix scan of max row values, i.e. remaining max weight in row
#endif
#if defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4)
        sum += rowval[i];
        sums[i] = sum;		// prefix sum
#endif
#if defined(DP3) || defined(DP4)
        hashsum[rowind[i]] = sums[i];	// hash sum
#endif
#if defined(DP7) || defined(DP8)
        hashsz[rowind[i]] = k++;
#endif
    }
#if defined(RS1) && !defined(RS3)
    if (rs1 < simT)
    {
#ifdef EXTRATIMES
        da_stopwctimer(params->timer_1);	// find candidates time
#endif
#ifndef NO_OUTPUT
        params->neighbors->rowptr[rid + 1] = params->neighbors->rowptr[rid];
#endif
        return;
    }
#endif
#endif

    // find minimum size
#ifdef SZ3
    minsize = eps / rwgts[rid];
#endif

    for (i = rowptr[rid + 1] - 1; i >= rowptr[rid]; i--)
    {
        k = rowind[i];
        myval = rowval[i];
        hashval[k] = myval;

#ifdef SZ3
        // remove those vectors from the inverted index that do
        // not have sufficient size.
        for (s = starts[k], j = starts[k]; j < ends[k]; j++)
        {
            cid = idxids[j];
            if (rsizes[cid]
#ifndef SZ1
                    * rwgts[cid]
#endif
                    >= minsize)
            break;
        }
        starts[k] += j - s;
#ifdef EXTRACOUNTS
        params->nPruneMinsize += j - s;
#endif
#endif

#ifdef RS1
#if defined(RS4)
        bnd = da_min(rs1, bsq);
#elif defined(RS2)
        bnd = da_min (rs1, rs2);
#else
        bnd = rs1;
#endif
#else
#if defined(RS4)
        bnd = bsq;
#elif defined(RS2)
        bnd = rs2;
#else
        bnd = 0;
        da_errexit (SIGERR,
                "Invalid pruning selected at compile time. One of RS1-RS4 must be selected.");
#endif
#endif

#if defined(IDXL2)
        rs2 -= (double) myval * myval;
        bsq = sqrt(rs2);
#elif defined(IDXLEN)
        rs2 -= (double) 0.5 *myval * myval;
#endif

#if defined(L2CV)
        hashlen[k] = lengths[i] = bsq;
#elif defined(LENCV)
        hashlen[k] = lengths[i] = rs2;
#endif


        for (j = starts[k]; j < ends[k]; j++) {
            cid = idxids[j];	// potential candidate from the inv index
            acc = accum[cid];

            if (acc > 0.0) {
                // accumulate
                accum[cid] = acc += myval * idxvals[j];
                // check l2-norm bound
#if defined(L2CG)
                if (acc + bsq * idxlens[j] < simTe) {
                    accum[cid] = -2;
#ifdef EXTRACOUNTS
                    params->nPruneLength++;
#endif
                }
#elif defined(LENCG)
                if (acc + rs2 + idxlens[j] < eps)
                {
                    accum[cid] = -2;
#ifdef EXTRACOUNTS
                    params->nPruneLength++;
#endif
                }
#endif
            } else if (bnd >= simT && acc > -2) {
                // accumulate
                if (acc == -1) {
                    accum[cid] = acc = myval * idxvals[j];
                    cands[ncands++] = cid;
                } else
                    accum[cid] = acc += myval * idxvals[j];
                // check l2-norm bound
#if defined(L2CG)
                if (acc + bsq * idxlens[j] < simTe) {
                    accum[cid] = -2;
#ifdef EXTRACOUNTS
                    params->nPruneLength++;
#endif
                }
#elif defined(LENCG)
                if (acc + rs2 + idxlens[j] < eps)
                {
                    accum[cid] = -2;
#ifdef EXTRACOUNTS
                    params->nPruneLength++;
#endif
                }
#endif
            }
        }
#ifdef RS1
#ifdef RS3
        rs1 -= myval * criwgts[k];
#else
        rs1 -= myval * cwgts[k];
#endif
#endif

    }
#ifdef EXTRATIMES
    da_stopwctimer(params->timer_1);	// find candidates time

    da_startwctimer(params->timer_2);	// process candidates time
#endif
    params->nCandidates += ncands;
    params->nDotProducts += l2apProcessCandidates(params, rid, docs, lengths,
            sums, endptr, ncands, cands, accum, hashval, hashwgt, hashlen,
            hashsum, hashsz, rwgts, rfiwgts, ps);
#ifdef EXTRATIMES
    da_stopwctimer(params->timer_2);	// process candidates time

    da_startwctimer(params->timer_1);	// find candidates time
#endif
    // reset hashval
    for (i = rowptr[rid]; i < endptr[rid]; i++)
        hashval[rowind[i]] = 0;

#ifdef EXTRATIMES
    da_stopwctimer(params->timer_1);	// find candidates time
#endif
}






idx_t
l2apProcessCandidates (params_t * params, idx_t rid, da_csr_t * docs,
        val_t * lengths, val_t * sums, ptr_t * endptr, idx_t ncands,
        idx_t * cands, accum_t * accum, val_t * hashval, val_t * hashwgt,
        val_t * hashlen, val_t * hashsum, idx_t * hashsz, val_t * rwgts,
        val_t * rfiwgts, val_t * ps)
{

    ssize_t i, j, k, p, r, ls, sz, cid, nneighbs = -1, dps = 0;
    val_t rwgt, rsum;
    double cval, simT, simTe;
    ptr_t *rowptr;
    idx_t *rowind;
    val_t *rowval;
    da_csr_t *neighbors;

    rowptr = docs->rowptr;	// index in rowind/rowval where each row starts
    rowind = docs->rowind;	// colid for each nnz
    rowval = docs->rowval;	// val for each nnz
#ifndef NO_OUTPUT
    if (params->nim)
    {
        neighbors = params->neighbors;
        // set current num neighbors for this row's search
        neighbors->rowptr[rid + 1] = neighbors->rowptr[rid];
        // grow neighborhood matrix if necessary
        if (neighbors->rowptr[rid] + ncands > params->nghnnz) {
            params->nghnnz += da_max(ncands, params->nghnnz / 2);
            da_csr_Grow(neighbors, params->nghnnz);
            params->nghinccnt++;
        }
    }
#endif

    simT = params->eps;
#if defined(L2CV)
    simTe = params->eps - LENBNDEPS;
#endif
#if defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7)
    rwgt = rwgts[rid];
#endif
#if defined(DP1) || defined(DP2)
    rsum = sums[rowptr[rid + 1] - 1];
#endif
#if defined(DP5) || defined(DP6)
    sz = rowptr[rid + 1] - rowptr[rid];
#endif

    for (i = 0; i < ncands; i++)
    {
        cid = cands[i];

        if (accum[cid] <= 0.0) {
            accum[cid] = -1;
            continue;
        }

        cval = accum[cid];
        accum[cid] = -1;

#ifdef PSCV
        if (cval + ((double) ps[cid]) < simT) {
#ifdef EXTRACOUNTS
            params->nPrunePscore++;
#endif
            continue;
        }
#endif

#ifdef DP1
        if (cval +
                ((double) da_min (rfiwgts[cid] * rsum, rwgt * sums[endptr[cid]])) <
                eps)
        {
#ifdef EXTRACOUNTS
            params->nPruneDotP++;
#endif
            continue;
        }
#endif

#ifdef DP5
        if (cval
                + ((double) da_min(endptr[cid] - rowptr[cid], sz) * rfiwgts[cid]
                        * rwgt) < simT) {
#ifdef EXTRACOUNTS
            params->nPruneDotP++;
#endif
            continue;
        }
#endif

        r = rowptr[cid];

#if defined(DP2) || defined(DP3) || defined(DP4) \
        || defined(DP6) || defined(DP7) || defined(DP8)

        // advance to next term in common
        for (j = endptr[cid] - 1, p = endptr[cid] - r; j >= r; j--, p--)
            if (hashval[rowind[j]] > 0.0)
                break;

        if (j < r) {
            params->nDotProducts1++; /** full dot product before verification **/
            goto check;
        }

#ifdef DP2
        if (cval +
                ((double)
                        da_min (rfiwgts[cid] * rsum, hashwgt[rowind[j]] * sums[j])) < eps)
        {
#ifdef EXTRACOUNTS
            params->nPruneDotP2++;
#endif
            continue;
        }
#elif defined(DP3)
        if (cval +
                ((double)
                        da_min (rfiwgts[cid] * hashsum[rowind[j]], rwgt * sums[j])) < eps)
        {
#ifdef EXTRACOUNTS
            params->nPruneDotP2++;
#endif
            continue;
        }
#elif defined(DP4)
        if (cval +
                ((double)
                        da_min (rfiwgts[cid] * hashsum[rowind[j]],
                                hashwgt[rowind[j]] * sums[j])) < eps)
        {
#ifdef EXTRACOUNTS
            params->nPruneDotP2++;
#endif
            continue;
        }
#elif defined(DP6)
        if (cval + ((double) da_min(p, sz) * rfiwgts[cid] * hashwgt[rowind[j]])
                < simT) {
#ifdef EXTRACOUNTS
            params->nPruneDotP2++;
#endif
            continue;
        }
#elif defined(DP7)
        if (cval +
                ((double) da_min (p, hashsz[rowind[j]]) * rfiwgts[cid] * rwgt) <
                eps)
        {
#ifdef EXTRACOUNTS
            params->nPruneDotP2++;
#endif
            continue;
        }
#elif defined(DP8)
        if (cval +
                ((double) da_min (p, hashsz[rowind[j]]) * rfiwgts[cid] *
                        hashwgt[rowind[j]]) < eps)
        {
#ifdef EXTRACOUNTS
            params->nPruneDotP2++;
#endif
            continue;
        }
#endif

#else
        j = endptr[cid] - 1;
#endif

#ifdef L2CV
        for (; j >= r; j--)
            if (hashval[rowind[j]] > 0.0) {
                cval += (double) hashval[rowind[j]] * rowval[j];
                if (j > r
                        && cval + ((double) hashlen[rowind[j]] * lengths[j])
                                < simTe) {
#ifdef EXTRACOUNTS
                    params->nPruneLength2++;
#endif
                    goto nexcid;
                }
            }
#elif defined(LENCV)
        for (; j >= r; j--)
        if (hashval[rowind[j]] > 0.0)
        {
            cval += (double) hashval[rowind[j]] * rowval[j];
            if (j > r
                    && cval + ((double) hashlen[rowind[j]]) +
                    ((double) lengths[j]) < eps)
            {
#ifdef EXTRACOUNTS
                params->nPruneLength2++;
#endif
                goto nexcid;
            }
        }
#else
        // partial dot product between rid and cid
        for (; j >= r; j--)
        if (hashval[rowind[j]] > 0.0)
        cval += hashval[rowind[j]] * rowval[j];
#endif
        dps++;

        check:			// check similarity value

        if (cval >= simT) {
            params->nSimPairs++;
            //add to neighbors matrix
            da_storeSim(params, rid, cid, cval);
        }

        nexcid: continue;

    }

    return dps;
}




/**
 * Index part of the vector after its search is complete
 */
void
l2apIndexRow (idx_t rid, params_t * params, da_csr_t * docs, ptr_t * endptr,
        da_invIdxJ_t * invIdx, val_t * rwgts, val_t * rfiwgts, val_t * cwgts,
        val_t * criwgts, val_t * hashwgt, val_t * ps)
{
    ssize_t j, k;
    val_t myval, maxrval;
    double b1, b2, b2l, sqb2, pscore, simT;
    ptr_t *rowptr = docs->rowptr;// index in rowind/rowval where each row starts
    idx_t *rowind = docs->rowind;	// colid for each nnz
    val_t *rowval = docs->rowval;	// val for each nnz
    simT = params->eps;
#ifdef PSCV
    pscore = 0.0;
#endif
#if defined(L2PS) || defined (IDXL2)
    b2 = sqb2 = b2l = 0.0;
#elif defined(LENPS) || defined (IDXLEN)
    b2 = 0.5;
#endif
#if ( defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7) ) && \
        !( defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8) )
    maxrval = 0.0;
#endif

    // index terms in row i
    myval = b1 = 0.0;
    for (endptr[rid] = rowptr[rid + 1], j = rowptr[rid]; j < endptr[rid]; j++)
    {
        myval = rowval[j];
        b1 += (double) myval * da_min(rwgts[rid], cwgts[rowind[j]]);
        // adjust b2 estimate of prefix similarity with all other vectors
#if defined(L2PS) || defined (IDXL2)
        b2 += (double) myval * myval;
        sqb2 = sqrt(b2);
#elif defined(LENPS) || defined (IDXLEN)
        b2 += (double) 0.5 *myval * myval;
#endif

        // check whether to start indexing
#ifdef L2PS
        if (da_min (b1, sqb2) >= simT)
            break;
#elif defined(LENPS)
        if (da_min (b1, b2) >= eps)
        break;
#else
        if (b1 >= eps)
        break;
#endif

#ifdef PSCV
        //update pscore
#ifdef L2PS
        pscore = da_min(b1, sqb2);
#elif defined(LENPS)
        pscore = da_min (b1, b2);
#else
        pscore = b1;
#endif
#endif

#if defined (L2PS)
        b2l = b2;
#endif

#if ( defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7) ) && \
        !( defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8) )
        // new f.i. max value
        if (myval > maxrval)
        maxrval = myval;
#endif

#if defined(RS1) && defined(RS3)
        // adjust inverted index max col weights
        if (myval > criwgts[rowind[j]])
        criwgts[rowind[j]] = myval;
#endif
    }
    // truncate the rest of the vector, since we're going to
    // index it.
    endptr[rid] = j;

    // update max forward index row val for row rid
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
    rfiwgts[rid] = j > rowptr[rid] ? hashwgt[rowind[j - 1]] : 0.0;
#elif defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7)
    rfiwgts[rid] = maxrval;
#endif

    //store the pscore for this row
#ifdef PSCV
    ps[rid] = pscore;
#endif

    // set b2 to the prefix norm up to current index but not including
#ifdef L2PS
    b2 = b2l;
    sqb2 = sqrt(b2);
#elif defined(LENPS)
    b2 -= 0.5;
#endif

    // index the remaining part of the vector
    for (; j < rowptr[rid + 1]; j++)
    {
        k = rowind[j];
#if defined (IDXL2)
        invIdx->lens[invIdx->ends[k]] = sqb2;
        b2 += (double) rowval[j] * rowval[j];
        sqb2 = sqrt(b2);
#elif defined(IDXLEN)
        invIdx->lens[invIdx->ends[k]] = b2;
        b2 += (double) 0.5 *rowval[j] * rowval[j];
#endif
        invIdx->ids[invIdx->ends[k]] = rid;
        invIdx->vals[invIdx->ends[k]] = rowval[j];
        invIdx->ends[k]++;
#if defined(RS1) && defined(RS3)
        if (rowval[j] > criwgts[k])
        criwgts[k] = rowval[j];
#endif
    }
}




/**
 * Reorder the document matrix in decreasing number of column nnz
 * and decreasing order of row max val. Also compute and store
 * necessary stats about the data
 * Store initial order in rperm/cperm.
 */
void l2apReorderDocs(params_t * params, da_csr_t ** docs, idx_t * rsizes,
        idx_t * csizes, val_t * rwgts, val_t * cwgts, idx_t * rperm,
        idx_t * cperm)
{
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
    for (j = 0; j < ncols; j++)
    {
        orderColNnz[j].key = j;
        orderColNnz[j].val = csizes[j];
    }
    da_iikvsortd(ncols, orderColNnz);	// sort columns
    for (j = 0; j < ncols; j++)
    {
        cperm[orderColNnz[j].key] = j;	// store column permutation
        csizes[j] = orderColNnz[j].val;	// new column size after permuting
    }

    // set new column ids in docs
    for (i = 0; i < nnz; i++)
        rowind[i] = cperm[rowind[i]];
    // get max weights for rows and columns + row sizes
    for (i = 0; i < nrows; i++)
    {
        rsizes[i] = (idx_t) (rowptr[i + 1] - rowptr[i]);	//record row nnzs
        for (j = rowptr[i]; j < rowptr[i + 1]; j++) {
            v = rowval[j];
            if (v > rwgts[i])
                rwgts[i] = v;	// store max val for row
            if (v > cwgts[rowind[j]])
                cwgts[rowind[j]] = v;	//store max val for col
        }
        //set up kv array to get new row order
        orderMaxVal[i].key = i;
        orderMaxVal[i].val = rwgts[i];
    }
    da_ivkvsortd(nrows, orderMaxVal);	// sort rows
    for (i = 0; i < nrows; i++)
    {
        rperm[i] = orderMaxVal[i].key;	// store new order
        rwgts[i] = orderMaxVal[i].val;	// new row weight after permuting
    }
    // permute rows in matrix
    for (nrowptr[0] = 0, nnz = 0, j = 0, i = 0; i < nrows; i++)
    {
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



/**
 * Reorder the document matrix in decreasing number of column nnzs.
 * Also compute and store necessary stats about the data.
 * Store initial order in cperm.
 */
void l2apReorderCols(da_csr_t * docs, idx_t * csizes, idx_t * cperm) {
    ssize_t i, j, k, nnz;
    idx_t nrows, ncols;
    da_iikv_t *orderColNnz = NULL;
    char pv = 0;
    ptr_t *rowptr;
    idx_t *rowind;
    val_t v, *rowval;

    nrows = docs->nrows;
    ncols = docs->ncols;
    rowptr = docs->rowptr;
    rowind = docs->rowind;
    nnz = rowptr[nrows];

    //record col nnzs
    da_iset(ncols, 0, csizes);
    for (i = 0; i < nnz; i++)
        csizes[rowind[i]]++;

    orderColNnz = da_iikvmalloc(ncols, "l2knnReorderCols: orderColNnz");
    //get new column order
    for (j = 0; j < ncols; j++)
    {
        orderColNnz[j].key = j;
        orderColNnz[j].val = csizes[j];
    }
    da_iikvsortd(ncols, orderColNnz);	// sort columns

    for (j = 0; j < ncols; j++)
    {
        if (cperm)
            cperm[orderColNnz[j].key] = j;	// store column permutation
        csizes[j] = orderColNnz[j].val;	// new column size after permuting
    }

    /* set new column ids in docs */
    for (i = 0; i < nnz; i++)
        rowind[i] = cperm[rowind[i]];

    da_csr_SortIndices(docs, DA_ROW);

    da_free((void **) &orderColNnz, LTERM);

}




/**
 * Reorder the document matrix in decreasing number of column nnzs.
 * Also compute and store necessary stats about the data.
 * Store initial order in cperm.
 */
void l2knnReorderCols(da_csr_t * docs, idx_t * csizes, idx_t * cperm) {
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
    cval = da_vsmalloc (ncols, 0.0, "l2knnReorderCols: cval");
#endif

    //record col nnzs
    da_iset(ncols, 0, csizes);
    for (i = 0; i < nnz; i++)
    {
#if COLSORTBY == COLSORTBY_MAX
        j = rowind[i];
        if (rowval[j] > cval[j])
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
    for (i = 0; i < ncols; ++i)
    if (csizes[i] > 0)
    cval[i] /= csizes[i];
#endif

#if COLSORTBY == COLSORTBY_MAX || COLSORTBY == COLSORTBY_AVG
    orderColVal = da_ivkvmalloc (ncols, "l2knnReorderCols: orderColNnz");
    //get new column order
    for (j = 0; j < ncols; j++)
    {
        orderColVal[j].key = j;
        orderColVal[j].val = cval[j];
    }
    da_ivkvsorti (ncols, orderColVal);	// sort columns

    for (j = 0; j < ncols; j++)
    {
        cperm[orderColVal[j].key] = j;	// store column permutation
        csizes[j] = orderColVal[j].val;// new column size after permuting
    }
#elif COLSORTBY == COLSORTBY_SIZE
    orderColNnz = da_iikvmalloc(ncols, "l2knnReorderCols: orderColNnz");
    //get new column order
    for (j = 0; j < ncols; j++)
    {
        orderColNnz[j].key = j;
        orderColNnz[j].val = csizes[j];
    }
    da_iikvsorti(ncols, orderColNnz);	// sort columns

    for (j = 0; j < ncols; j++)
    {
        cperm[orderColNnz[j].key] = j;	// store column permutation
        csizes[j] = orderColNnz[j].val;	// new column size after permuting
    }
#elif COLSORTBY == COLSORTBY_RAND
    da_irandArrayPermute (ncols, cperm, 10 * ncols, 1);
    tsizes = da_imalloc (ncols, "l2knnReorderCols: tsizes");
    for (j = 0; j < ncols; j++)
    {
        tsizes[j] = csizes[cperm[j]];	// new column size after permuting
    }
    da_icopy (ncols, tsizes, csizes);
#else
    da_errexit ("Invalid COLSORTBY criteria in defs.h.");
#endif

    /* set new column ids in docs */
    for (i = 0; i < nnz; i++)
        rowind[i] = cperm[rowind[i]];

    da_csr_SortIndices(docs, DA_ROW);

    da_free((void **) &orderColNnz, &orderColVal, &cval, &tsizes, LTERM);

}
