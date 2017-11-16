/*!
 \file  l2ap-c.c
 \brief This file contains related functions for the Tanimoto version of l2ap

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
void l2apFindMatchesTan(
        idx_t const rid,
        idx_t const minid,
        val_t const maxnrm,
        params_t * const params,
        da_csr_t const * const fidx,
        da_csr_t const * const iidx,
        da_csr_t const * const docs,
        idx_t const * const rsizes,
        val_t const * const rwgts,
        val_t const * const rfiwgts,
        val_t const * const cwgts,
        val_t const * const ps,
        val_t const * const plen,
        char * const mark,
        idx_t * const cands,
        accum_t * const accum,
        val_t * const hashval,
        val_t * const hashlen );
idx_t l2apProcessCandidatesTan(
        idx_t const rid,
        idx_t const ncands,
        params_t * const params,
        da_csr_t const * const fidx,
        val_t const * const rnorms,
        idx_t const * const rsizes,
        val_t const * const rwgts,
        val_t const * const rfiwgts,
        val_t const * const ps,
        char * const mark,
        idx_t * const cands,
        accum_t * const accum,
        val_t const * const hashval,
        val_t const * const hashlen );
void l2apTanCreateIndices(
        params_t const * const params,
        da_csr_t * docs,
        da_csr_t * iidx,
        da_csr_t * fidx,
        val_t * const rfiwgts,
        val_t * const cwgts,
        val_t * const ps,
        ptr_t * const eptr);
void l2apIndexRowTan(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdxJ_t *invIdx,
        val_t *rwgts, val_t *rfiwgts, val_t *cwgts, val_t *hashwgt, val_t *ps);
void l2apReorderDocsTan(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
        val_t *cwgts, idx_t *rperm, idx_t *cperm);

/**
 * Main entry point to l2ap for Tanimoto.
 */
void l2apFindNeighborsTan(params_t *params){

    ssize_t i, j;
    idx_t minid, progressInd;
    double amin, lrnrm, v;
    da_csr_t *docs;

    if(params->norm > 0){
        printf("Warning: setting -norm=0 since -sim=tan.\n");
        params->norm = 0;
    }

    // allocate memory
    simSearchSetup(params);

    params->eps = (2.0*(double)params->simT)/(1.0+(double)params->simT) - TANEPS/2;
    double const b = params->b = 1.0/((double)params->eps);
    params->b2 = b*b;
    double const a = params->a = b + sqrt(b*b-1.0);         // alpha/2

    docs   = params->docs;
    idx_t const nrows  = docs->nrows;   // num rows
    idx_t const ncols  = docs->ncols;   // num cols
    progressInd = ceil(nrows/(float)10);

    idx_t * const rsizes  = da_imalloc(nrows, "l2apFindNeighborsTan: rsizes");
    idx_t * const csizes  = da_ismalloc(ncols, 0.0, "l2apFindNeighborsTan: csizes");
    idx_t * const rperm   = da_imalloc(nrows, "l2apFindNeighborsTan: rperm");
    idx_t * const cperm   = da_imalloc(ncols, "l2apFindNeighborsTan: cperm");
    val_t * const rwgts   = da_vsmalloc(nrows, 0.0, "l2apFindNeighborsTan: rwgts");
    val_t * const rfiwgts = da_vmalloc(nrows, "l2apFindNeighborsTan: rfiwgts");
    val_t * const cwgts   = da_vsmalloc(ncols, 0.0, "l2apFindNeighborsTan: cwgts");
    val_t * const hashval = da_vsmalloc(ncols, 0, "l2apFindNeighborsTan: hashval");      // hash values of x
#if defined(L2CV) || defined(TL2)
    val_t * const hashlen = da_vsmalloc(ncols, 0, "l2apFindNeighborsTan: hashlen");      // hash suffix lengths of x
#endif
#ifdef PSCV
    val_t * const ps      = da_vmalloc(nrows, "l2ap2FindNeighbors: ps");        // score prior to indexing threshold for each row
#endif

    // reorder matrix
    l2apReorderDocsTan(params, &docs, rsizes, csizes, rwgts, cwgts, rperm, cperm);

    // compute prefix l2-norms in original space and then normalize the rows
#ifdef TRS
    val_t * const plen    = da_vsmalloc(docs->rowptr[nrows], 0.0, "l2apFindNeighborsTan: cwgts");
    for(i=0; i < nrows; i++){
        for(v=0.0, j=docs->rowptr[i]; j < docs->rowptr[i+1]; ++j){
            v += docs->rowval[j] * docs->rowval[j];
            plen[j] = sqrt(v);
        }
    }
#else
    val_t * const plen    = NULL;
#endif
    da_csr_Normalize(docs, DA_ROW, 2);


    // create indices
    da_csr_t * const fidx = da_csr_Create();
    da_csr_t * const iidx = da_csr_Create();
    ptr_t * const eptr = da_psmalloc(nrows+1, 0, "l2apFindNeighborsTan: eptr");
    l2apTanCreateIndices(params, docs, iidx, fidx, rfiwgts, cwgts, ps, eptr);

    ptr_t const * const ptr = docs->rowptr; // index in rowind/rowval where each row starts
    val_t const * const rnorms = docs->rnorms;
    minid  = 0;

    idx_t * const cands  = da_imalloc(nrows, "l2apFindNeighborsTan: candidates");
    char * const mark = da_csmalloc(nrows, S_INIT, "l2apFindNeighborsTan: invIdx->mark");
    accum_t * const accum = gk_malloc(nrows*sizeof(accum_t), "l2apFindNeighborsTan: invIdx->accum");
    for(i=0; i < nrows; i++){
        accum[i] = 0.0;
    }

#ifdef EXTRACOUNTS
    // find the average length of rows & columns
    idx_t const mrlen = ceil(da_isum(nrows, rsizes, 1) / (val_t) nrows);
    idx_t const mclen = ceil(da_isum(ncols, csizes, 1) / (val_t) ncols);
    if(params->verbosity > 0)
        printf("mrlen: %d, mclen: %d.\n", mrlen, mclen);
#endif

    if(params->verbosity > 0)
        printCompileChoices();

    // for each x \in V do
    gk_startwctimer(params->timer_3); // find neighbors time
    if(params->verbosity > 0)
        printf("Progress Indicator (no. of records): ");

    lrnrm = rnorms[0];
    for(i=0; i < nrows; i++){

        // find mind and max id
        val_t const rnrm = rnorms[i];

        for(j=minid; j < i; ++j){
            if(rnorms[j] * a >= rnrm){
                minid = j;
                break;
            }
        }

        // find matches for row i
        l2apFindMatchesTan(i, minid, lrnrm, params, fidx, iidx, docs, rsizes, rwgts, rfiwgts, cwgts,
                ps, plen, mark, cands, accum, hashval, hashlen);

        lrnrm = rnrm;

        if ( params->verbosity > 0 && i % progressInd == 0 ) {
            printf("%d%%..", (int) (100*i/(float)nrows));
            fflush(stdout);
        }
    }

    if(params->verbosity > 0)
        printf("100%%\n");


#ifdef EXTRACOUNTS
    params->indexSize = iidx->rowptr[iidx->nrows];
    params->nCandidates -= params->nPruneTanLen;
#endif

    gk_stopwctimer(params->timer_3); // find neighbors time

    simSearchFinalize(params);
    da_csr_FreeAll((da_csr_t**)&iidx, &fidx, LTERM);
    gk_free((void**)&rsizes, &csizes, &rwgts, &rfiwgts, &cwgts, &eptr,
            &cands, &hashval, &hashlen, &ps, &mark, &accum, &plen, LTERM);

}



/**
 * Find matches for row id rid
 */
void l2apFindMatchesTan(
        idx_t const rid,
        idx_t const minid,
        val_t const maxnrm,
        params_t * const params,
        da_csr_t const * const fidx,
        da_csr_t const * const iidx,
        da_csr_t const * const docs,
        idx_t const * const rsizes,
        val_t const * const rwgts,
        val_t const * const rfiwgts,
        val_t const * const cwgts,
        val_t const * const ps,
        val_t const * const plen,
        char * const mark,
        idx_t * const cands,
        accum_t * const accum,
        val_t * const hashval,
        val_t * const hashlen )
{

#ifdef EXTRATIMES
    gk_startwctimer(params->timer_1); // find candidates time
#endif
    ptr_t i, j;
    idx_t ncands;
    accum_t simT, cval, bnd;

    ptr_t const * const ptr = docs->rowptr;
    idx_t const * const ind = docs->rowind;
    val_t const * const val = docs->rowval;
    val_t const * const len = docs->rvols;

    ptr_t const * const iptr = iidx->rowptr;
    idx_t const * const iind = iidx->rowind;
    val_t const * const ival = iidx->rowval;
    val_t const * const ilen = iidx->rvols;

    ncands  = 0;            // number of candidates for this doc
    bnd = 1.0;

    simT  = params->eps;

#ifdef TRS
    const val_t rnrm = docs->rnorms[rid];
    const val_t minnrm = docs->rnorms[0];
    const val_t b2 = 0.5 * params->eps * (rnrm * rnrm + minnrm*minnrm)/maxnrm;
#endif

    for ( i=ptr[rid+1]-1; i >= ptr[rid]; i-- )
    {
        if(bnd < simT){
            #ifdef EXTRACOUNTS
            params->nCountRS++;
            #endif
            break;
        }
#ifdef TRS
        else
        if(bnd * plen[i] < b2){
            #ifdef EXTRACOUNTS
            params->nCountTanRS++;
            #endif
            break;
        }
#endif

        val_t const l = len[i];
        idx_t const c = ind[i];
        val_t const v = val[i];
        hashval[c] = v;
        hashlen[c] = bnd = l;


        for ( j = iptr[c]; j < iptr[c+1]; j++ ){
            idx_t const cid = iind[j];  // potential candidate from the inv index
            if(cid >= rid){
                break;
            }
            if ( mark[cid] == S_START ){
                // accumulate
                accum[cid] += v * ival[j];
                // check l2-norm bound
                #if defined(L2CG)
                if(accum[cid] + l * ilen[j] < simT){
                    mark[cid] = S_PRUNED;
                    #ifdef EXTRACOUNTS
                    params->nPruneLength++;
                    #endif
                }
                #endif

            } else if ( mark[cid] == S_INIT){
                if(cid < minid){
                    // prune based on Tanimoto bound on orig. vector length
                    #ifdef EXTRACOUNTS
                    mark[cid] = S_PRUNED;
                    params->nPruneTanLen++;
                    cands[ncands++] = cid;
                    #endif
                    continue;
                }

                // accumulate
                accum[cid] = v * ival[j];
                cands[ncands++] = cid;
                mark[cid] = S_START;

                // check l2-norm bound
                #if defined(L2CG)
                if(accum[cid] + l * ilen[j] < simT){
                    mark[cid] = S_PRUNED;
                    #ifdef EXTRACOUNTS
                    params->nPruneLength++;
                    #endif
                }
                #endif
            }
        }

    }



    for ( ; i >= ptr[rid]; i-- )
    {
        val_t const l = len[i];
        idx_t const c = ind[i];
        val_t const v = val[i];
        hashval[c] = v;
        hashlen[c] = l;

        for ( j = iptr[c]; j < iptr[c+1]; j++ ){
            idx_t const cid = iind[j];  // potential candidate from the inv index
            if(cid >= rid){
                break;
            }
            if ( mark[cid] == S_START ){
                // accumulate
                accum[cid] += v * ival[j];
                // check l2-norm bound
                #if defined(L2CG)
                if(accum[cid] + l * ilen[j] < simT){
                    mark[cid] = S_PRUNED;
                    #ifdef EXTRACOUNTS
                    params->nPruneLength++;
                    #endif
                }
                #endif

            }
        }

    }



#ifdef EXTRATIMES
    gk_stopwctimer(params->timer_1); // find candidates time

    gk_startwctimer(params->timer_2); // process candidates time
#endif
    params->nCandidates  += ncands;
    params->nDotProducts += l2apProcessCandidatesTan(rid, ncands, params, fidx, docs->rnorms, rsizes,
            rwgts, rfiwgts, ps, mark, cands, accum, hashval, hashlen);
#ifdef EXTRATIMES
    gk_stopwctimer(params->timer_2); // process candidates time

    gk_startwctimer(params->timer_1); // find candidates time
#endif
    // reset hashval
    for ( i=ptr[rid]; i < ptr[rid+1]; i++ ){
        hashval[ind[i]] = 0;
    }

#ifdef EXTRATIMES
    gk_stopwctimer(params->timer_1); // find candidates time
#endif
}



idx_t l2apProcessCandidatesTan(
        idx_t const rid,
        idx_t const ncands,
        params_t * const params,
        da_csr_t const * const fidx,
        val_t const * const rnorms,
        idx_t const * const rsizes,
        val_t const * const rwgts,
        val_t const * const rfiwgts,
        val_t const * const ps,
        char * const mark,
        idx_t * const cands,
        accum_t * const accum,
        val_t const * const hashval,
        val_t const * const hashlen )
{

    idx_t i;
    ptr_t j, p, r, dps = 0;
    val_t ubnd;
    accum_t cval;
    da_csr_t *neighbors;

    ptr_t const * const fptr = fidx->rowptr;
    idx_t const * const find = fidx->rowind;
    val_t const * const fval = fidx->rowval;
    val_t const * const flen = fidx->rvols;


    val_t const rnrm   = rnorms[rid];


    accum_t const simT = params->eps;
    #if defined(TL1) || defined (TL2)
    accum_t const b = params->b;
    accum_t const b2 = params->b2;
    #endif

    if(params->nim){
        neighbors = params->neighbors;
        // set current num neighbors for this row's search
        neighbors->rowptr[rid+1] = neighbors->rowptr[rid];
        // grow neighborhood matrix if necessary
        if(neighbors->rowptr[rid] + ncands > params->nghnnz){
            params->nghnnz += gk_max(ncands, params->nghnnz/2);
            da_csr_Grow(neighbors, params->nghnnz);
            params->nghinccnt++;
        }
    }

#if defined(L2CV)
    accum_t const simTe = simT-LENBNDEPS;
#endif
#if defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7)
    val_t const rwgt  = rwgts[rid];
#endif
#if defined(DP5) || defined(DP6)
    val_t const sz    = rsizes[rid];
#endif

    for ( i=0; i < ncands; i++ )
    {
        idx_t const cid = cands[i];

        if(mark[cid] == S_PRUNED){
            mark[cid] = S_INIT;
            continue;
        }

        cval = accum[cid];
        accum[cid] = 0.0;
        mark[cid] = S_INIT;


#ifdef PSCV
        if ( (ubnd = cval + ps[cid]) < simT ){
            #ifdef EXTRACOUNTS
                params->nPrunePscore++;
            #endif
            continue;
        }

    #ifdef TL1
        if(rnorms[cid] * ( ubnd*b + sqrt(ubnd*ubnd*b2 - 1.0) ) < rnrm){
            #ifdef EXTRACOUNTS
                params->nPruneTanLenCV++;
            #endif
            continue;
        }
    #endif
#endif

#ifdef DP5
        if( cval + gk_min(fptr[cid+1]-fptr[cid], sz) * rfiwgts[cid] * rwgt < simT){
            #ifdef EXTRACOUNTS
                params->nPruneDotP++;
            #endif
            continue;
        }
#endif



#if defined(TL2)

        // advance to next term in common
        r = fptr[cid];
        for(j=fptr[cid+1]-1, p=fptr[cid+1]-1-r; j >= r; j--, p--){
            if(hashval[find[j]] > 0.0){
                break;
            }
        }

        if(j < r)
            goto check;

    #ifdef TL2
        if(hashval[find[j]] > 0.0){
            ubnd = cval + hashval[find[j]] * fval[j] + hashlen[find[j]] * flen[j];

            if(rnorms[cid] * ( ubnd*b + sqrt(ubnd*ubnd*b2 - 1.0) ) < rnrm){
                #ifdef EXTRACOUNTS
                    params->nPruneTanLenCV++;
                #endif
                continue;
            }
        }
    #endif

#else
        r = fptr[cid];
        j = fptr[cid+1]-1;
#endif

#ifdef L2CV
        for( ; j >= r; j--){
            if(hashval[find[j]] > 0.0){
                cval += hashval[find[j]] * fval[j];
                if(j > r && cval + hashlen[find[j]] * flen[j] < simTe){
                    #ifdef EXTRACOUNTS
                        params->nPruneLength2++;
                    #endif
                    goto nexcid;
                }
            }
        }
#else
        // partial dot product between rid and cid
        for( ; j >= r; j--)
            if(hashval[find[j]] > 0.0)
                cval += hashval[find[j]] * fval[j];
#endif

        check:  // check similarity value

        dps++;

        //convert to tanimoto from cosine and check initial threshold
        cval /= ( (rnrm*rnrm +
                   rnorms[cid]*rnorms[cid])/(rnrm*rnorms[cid]) - cval);
        if(cval >= params->simT - TANEPS){
            //add to neighbors matrix
            params->nSimPairs++;
            da_storeSim(params, rid, cid, cval);
        }

        nexcid:
            continue;

    }

    return dps;
}


/**
 * Index part of the vector after its search is complete
 */
void l2apIndexRowTan(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdxJ_t *invIdx,
        val_t *rwgts, val_t *rfiwgts, val_t *cwgts, val_t *hashwgt, val_t *ps){
    ssize_t j, k;
    val_t myval, maxrval, t;
    accum_t b1, b2, b2l, sqb2, pscore, simT;
    ptr_t *rowptr = docs->rowptr; // index in rowind/rowval where each row starts
    idx_t *rowind = docs->rowind; // colid for each nnz
    val_t *rowval = docs->rowval; // val for each nnz


    simT = params->eps;


#ifdef PSCV
    pscore = 0.0;
#endif
#if defined(L2PS) || defined (IDXL2)
    b2 = sqb2 = b2l = 0.0;
#endif
#if ( defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7) ) && \
    !( defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8) )
    maxrval = 0.0;
#endif

    // index terms in row i
    myval = b1 = 0.0;
    for ( endptr[rid] = rowptr[rid+1], j = rowptr[rid]; j < endptr[rid]; j++ ) {
        myval = rowval[j];
        b1 += myval * cwgts[rowind[j]];
        // adjust b2 estimate of prefix similarity with all other vectors
#if defined(L2PS) || defined (IDXL2)
        b2 += myval * myval;
        sqb2 = sqrt(b2);
#endif

        // check whether to start indexing
#ifdef L2PS
        if(gk_min(b1,sqb2) >= simT)
            break;
#else
        if(b1 >= simT )
            break;
#endif

#ifdef PSCV
        //update pscore
    #ifdef L2PS
        pscore = gk_min(b1,sqb2);
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
        if( myval > maxrval )
            maxrval = myval;
#endif
    }
    // truncate the rest of the vector, since we're going to
    // index it.
    endptr[rid] = j;

    // update max forward index row val for row rid
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
    rfiwgts[rid] = j > rowptr[rid] ? hashwgt[rowind[j-1]] : 0.0;
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
#endif

    // index the remaining part of the vector
    for ( ; j < rowptr[rid+1]; j++ ){
        k = rowind[j];
#if defined (IDXL2)
        invIdx->lens[ invIdx->ends[k] ] = sqb2;
        b2 += rowval[j] * rowval[j];
        sqb2 = sqrt(b2);
#endif
        invIdx->ids[ invIdx->ends[k] ] = rid;
        invIdx->vals[ invIdx->ends[k] ] = rowval[j];
        invIdx->ends[k]++;
    }
}


/**
 * Create the inverted and forward indices. iidx and fidx are assumed empty csr structures
 */
void l2apTanCreateIndices(
        params_t const * const params,
        da_csr_t * docs,
        da_csr_t * iidx,
        da_csr_t * fidx,
        val_t * const rfiwgts,
        val_t * const cwgts,
        val_t * const ps,
        ptr_t * const eptr)
{
    idx_t i, k;
    ptr_t j, nnz, nnz2;
    val_t pscore, max, v;
    accum_t b, b2, sqb2;

    val_t const eps = params->eps;
    idx_t const nrows = docs->nrows;
    idx_t const ncols = docs->ncols;
    ptr_t const * const ptr = docs->rowptr;
    idx_t const * const ind = docs->rowind;
    val_t const * const val = docs->rowval;
    val_t * const len = docs->rvols = da_vmalloc(ptr[nrows], "l2apTanCreateIndices: len"); /* prefix l2 norms of normalized vectors */

    /* find the indexing point for each item and compute suffix lengths */
    for(nnz=0, i=0; i < nrows; ++i){
        pscore = max = b = b2 = sqb2 = 0.0;
        eptr[i] = ptr[i+1];
        for(j=ptr[i]; j < ptr[i+1]; ++j){
            v = val[j];

            // compute suffix lengths
            len[j] = sqb2;

            b += v * cwgts[ind[j]];
            b2 += v * v;
            sqb2 = sqrt(b2);

            if(gk_min(b, sqb2) >= eps){
                break;
            }

            pscore = gk_min(b, sqb2);

            if(v > max){
                max = v;
            }
        }
        eptr[i] = j;
        rfiwgts[i] = max;
        ps[i] = pscore;
        nnz += j - ptr[i];

        // compute remaining suffix lengths
        for( j++ ; j < ptr[i+1]; ++j){
            len[j] = b2 < 1.0 ? sqrt(b2) : 1.0;
            b2 += val[j] * val[j];
        }
    }

    /* allocate space */
    fidx->nrows = nrows;
    fidx->ncols = ncols;
    ptr_t * const fptr = fidx->rowptr = da_psmalloc(nrows+1, 0, "l2apTanCreateIndices: fptr");
    idx_t * const find = fidx->rowind = da_imalloc(nnz, "l2apTanCreateIndices: find");
    val_t * const fval = fidx->rowval = da_vmalloc(nnz, "l2apTanCreateIndices: fval");
    val_t * const flen = fidx->rvols  = da_vmalloc(nnz, "l2apTanCreateIndices: flen");

    nnz2 = ptr[nrows] - nnz;
    iidx->nrows = ncols;
    iidx->ncols = nrows;
    ptr_t * const iptr = iidx->rowptr = da_psmalloc(ncols+1, 0, "l2apTanCreateIndices: iptr");
    idx_t * const iind = iidx->rowind = da_imalloc(nnz2, "l2apTanCreateIndices: iind");
    val_t * const ival = iidx->rowval = da_vmalloc(nnz2, "l2apTanCreateIndices: ival");
    val_t * const ilen = iidx->rvols  = da_vmalloc(nnz2, "l2apTanCreateIndices: ilen");

    /* create indices */
    for (i=0; i<nrows; i++) {
        fptr[i] = eptr[i]-ptr[i];
        for (j=eptr[i]; j < ptr[i+1]; j++){
            iptr[ind[j]]++;
        }
    }
    MAKECSR(i, nrows, fptr);
    ASSERT(fptr[nrows] == nnz);
    MAKECSR(i, ncols, iptr);
    ASSERT(iptr[ncols] == nnz2);
    for (i=0; i<nrows; i++) {
        k = eptr[i]-ptr[i];
        if(k > 0){
            da_icopy(k, docs->rowind+ptr[i], find+fptr[i]);
            da_vcopy(k, docs->rowval+ptr[i], fval+fptr[i]);
            da_vcopy(k, len+ptr[i], flen+fptr[i]);
        }
        for (j=eptr[i] ; j < ptr[i+1]; j++){
            k = ind[j];
            iind[iptr[k]]   = i;
            ival[iptr[k]]   = val[j];
            ilen[iptr[k]++] = len[j];
        }
    }
    SHIFTCSR(i, ncols, iptr);

}




/**
 * Reorder the document matrix in decreasing number of column nnz
 * and increasing row length order. Also compute and store
 * necessary stats about the data
 * Store initial order in rperm/cperm.
 */
void l2apReorderDocsTan(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
        val_t *cwgts, idx_t *rperm, idx_t *cperm)
{
    ssize_t i, j, k, nnz;
    idx_t nrows, ncols;
    da_iikv_t *orderColNnz = NULL;
    da_ivkv_t *orderL2Norm = NULL;
    char pv = 0;
    ptr_t *rowptr, *nrowptr;
    idx_t *rowind, *nrowind;
    val_t v, *rowval, *nrowval, *rnorms, *trwgts=NULL;

    nrows  = (*docs)->nrows;
    ncols  = (*docs)->ncols;
    rowptr = (*docs)->rowptr;
    rowind = (*docs)->rowind;
    rowval = (*docs)->rowval;
    nnz    = rowptr[nrows];

    // for tan similarity, first compute and store norms, than normalize
    da_csr_ComputeSquaredNorms(*docs, DA_ROW);
    rnorms = (*docs)->rnorms;
    // now normalize the rows
    for(i=0; i < nrows; ++i){
        rnorms[i] = sqrt(rnorms[i]);
    }
    trwgts = da_vmalloc(nrows, "l2apReorderDocsTan: rwgts");

    //record col nnzs
    da_iset(ncols, 0, csizes);
    for (i=0; i<nnz; i++)
        csizes[rowind[i]]++;

    //assign memory for ordering
    orderL2Norm = da_ivkvmalloc(nrows, "l2apReorderDocsTan: orderMaxVal");
    nrowptr = da_pmalloc(nrows+1, "l2apReorderDocsTan: nrowptr");
    nrowind = da_imalloc(nnz, "l2apReorderDocsTan: nrowind");
    nrowval = da_vmalloc(nnz, "l2apReorderDocsTan: nrowval");
    orderColNnz = da_iikvmalloc(ncols, "l2apReorderDocsTan: orderColNnz");

    //get new column order
    for(j=0; j < ncols; j++){
        orderColNnz[j].key = j;
        orderColNnz[j].val = csizes[j];
    }
    if(params->mode == MODE_AP2)
        da_iikvsorti(ncols, orderColNnz); // sort columns
    else
        da_iikvsortd(ncols, orderColNnz); // sort columns
    for(j=0; j < ncols; j++){
        cperm[orderColNnz[j].key] = j;  // store column permutation
        csizes[j] = orderColNnz[j].val; // new column size after permuting
    }

    // set new column ids in docs
    for (i=0; i<nnz; i++)
        rowind[i] = cperm[rowind[i]];
    // get max weights for rows and columns + row sizes
    for (i=0; i<nrows; i++) {
        rsizes[i] = (idx_t) (rowptr[i+1] - rowptr[i]); //record row nnzs
        for (j=rowptr[i]; j<rowptr[i+1]; j++) {
            v = rowval[j];
            if(v > rwgts[i])
                rwgts[i] = v; // store max val for row
            if(v > cwgts[rowind[j]])
                cwgts[rowind[j]] = v; //store max val for col
        }
        //set up kv array to get new row order
        orderL2Norm[i].key = i;
        orderL2Norm[i].val = rnorms[i];
    }
    da_ivkvsorti(nrows, orderL2Norm); // sort rows
    for(i=0; i < nrows; i++){
        rperm[i] = orderL2Norm[i].key; // store new order
        rnorms[i] = orderL2Norm[i].val; // new row weight after permuting
    }
    // permute rows in matrix
    for (nrowptr[0] = 0, nnz=0, j=0, i=0; i<nrows; i++) {
        k = rperm[i];
        da_icopy(rowptr[k+1]-rowptr[k], rowind+rowptr[k], nrowind+nnz);
        da_vcopy(rowptr[k+1]-rowptr[k], rowval+rowptr[k], nrowval+nnz);
        rsizes[i] = rowptr[k+1]-rowptr[k];
        nnz += rsizes[i];
        nrowptr[++j] = nnz;
    }
    gk_free((void**)&(*docs)->rowptr, &(*docs)->rowind, &(*docs)->rowval, LTERM);
    // permute norms
    for(i=0; i < nrows; ++i){
        trwgts[i] = rwgts[rperm[i]];
    }
    da_vcopy(nrows, trwgts, rwgts);
    (*docs)->rowptr = nrowptr;
    (*docs)->rowind = nrowind;
    (*docs)->rowval = nrowval;

    da_csr_SortIndices((*docs), DA_ROW);

    gk_free((void**)&orderColNnz, &orderL2Norm, &cperm, &trwgts, LTERM);

    params->rperm = rperm;

}

