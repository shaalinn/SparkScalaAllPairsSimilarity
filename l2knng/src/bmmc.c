/*!
 \file  bmmc.c
 \brief This file contains Block-Max Maxscore (BMM) related functions __WITH block compression__.
         BMM is adapted to work with the cosine function, following the method in the paper below.
         Then, a k-nn query is executed with BMM for each document in the dataset to obtain the final
         KNNG.

@inproceedings{Dimopoulos:2013:OTD:2433396.2433412,
 author = {Dimopoulos, Constantinos and Nepomnyachiy, Sergey and Suel, Torsten},
 title = {Optimizing Top-k Document Retrieval Strategies for Block-max Indexes},
 booktitle = {Proceedings of the Sixth ACM International Conference on Web Search and Data Mining},
 series = {WSDM '13},
 year = {2013},
 isbn = {978-1-4503-1869-3},
 location = {Rome, Italy},
 pages = {113--122},
 numpages = {10},
 url = {http://doi.acm.org/10.1145/2433396.2433412},
 doi = {10.1145/2433396.2433412},
 acmid = {2433412},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {block-max inverted index, docid-oriented block-max index, early termination, top-k query processing},
}

 \author David C. Anastasiu
 */

#include "includes.h"
#include "pfor.h"

#define COMPR_BS   64  // block size for compression

// some helpful macros  //

/* what block is document in - doc block*/
#define DBL(did, bxp)\
    (did >> bxp)

/* map document to some aggregate block data */
#define MAP(arr, did, bxp)\
    (arr[(did >> bxp)])

/* what is the first potential doc in next block - next block doc */
#define NEXTBLD(did, bxp)\
    (((did >> bxp) + 1) << bxp)


// forward declarations
inline float getMaxScore(da_bil_t *eptr, uint did);
inline float getMaxDid(da_bil_t *eptr, uint did);
inline uint getDocId(da_bil_t *eptr);
inline uint getFreq(da_bil_t *eptr);
inline uint decompressBlock(da_bil_t *eptr);
inline uint decompressBlockNum(da_bil_t *eptr, uint bid);
inline uint goToDocId(da_bil_t *eptr, uint did);
inline uint nextLiveBlock(da_bil_t **eptrs, idx_t sz, uint did);
inline uint nextDocId(da_bil_t **eptrs, idx_t sz);
inline void listsToDocId(da_bil_t **eptrs, idx_t sz, uint did);
inline uint nextDoc(da_bil_t *eptr);
uint* compressList(uint *where, uint *ids, uint *freq, const uint sz, const uint prefix);
uint decompressList(uint *what, uint *ids, uint *freq);
inline float computeSim(idx_t rid, idx_t did, val_t rnorm, val_t cnorm, idx_t sz, da_bil_t **eptrs);
void da_bmmcGetSimilarRows(uint rid, da_csr_t *docs, da_iapq_t *knn, da_bil_t **idx, uint *idz,
        da_bil_t **qptrs, val_t *pmscore, val_t *bmscore, params_t *params);
da_bil_t** createBlockIndexC(da_csr_t *docs, val_t *idfs, val_t *norms, uint *ids, uint *freqs);
void freeBlockIndexC(da_bil_t ***idxp, idx_t sz);

/** block sizes for max doc id in list -- e.g., if max 2^16 <= docId < 2^17, block size for list is BLSIZES[17] = 8. */
static char BLSIZES[] = {
        10,   /* 2^0 */
        10,   /* 2^1 */
        10,   /* 2^2 */
        10,   /* 2^3 */
        10,   /* 2^4 */
        10,   /* 2^5 */
        10,   /* 2^6 */
        10,   /* 2^7 */
        10,   /* 2^8 */
        10,   /* 2^9 */
        10,   /* 2^10 */
        6,    /* 2^11 */
        6,    /* 2^12 */
        7,    /* 2^13 */
        7,    /* 2^14 */
        7,    /* 2^15 */
        8,    /* 2^16 */
        8,    /* 2^17 */
        7,    /* 2^18 */
        6,    /* 2^19 */
        6,    /* 2^20 */
        6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};

#define MAXBLSZ (1<<10)

/**
 * Main entry point to BMMC.
 */
void bmmcFindNeighbors(params_t *params)
{

    ssize_t i, j, k, nneighbs;
    size_t rid, nnz;
    idx_t nrows, ncols, progressInd, pct, mrlen;
    ptr_t *rowptr, *colptr, *nptr;
    idx_t *rowind, *colind, *nind, *cfreq;
    val_t *rowval, *colval, *nval, *mscore, *idfs, *norms, *pmscore, *bmscore;
    da_csr_t *docs, *neighbors=NULL;
    da_bil_t **qptrs = NULL;
    da_iapq_t *knn = NULL;
    da_bil_t **idx = NULL;
    uint *ids, *freqs;
    double norm;
    float idf, sc, msc;

    if(params->fpout && params->fmtWrite != DA_FMT_IJV)
        params->nim = 1;

    docs    = params->docs;
    nrows   = docs->nrows;  // num rows
    rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
    rowind  = docs->rowind; // colid for each nnz
    rowval  = docs->rowval; // val for each nnz
    norms   = docs->rnorms;
    nnz     = rowptr[nrows];

    da_progress_init_steps(pct, progressInd, nrows, 100);

    da_startwctimer(params->timer_3); // overall knn graph construction time

    da_startwctimer(params->timer_7); // indexing time

    // compact the column space
    da_csr_SortIndices(docs, DA_ROW);
    da_csr_CompactColumns(docs);

    // create inverted index and sort inv lists in doc id order
    if(!docs->colptr)
        da_csr_CreateIndex(docs, DA_COL);
    da_csr_SortIndices(docs, DA_COL);

    ncols   = docs->ncols;  // num cols
    colptr  = docs->colptr; // index in colind/colval/colids where each column starts
    colind  = docs->colind;
    colval  = docs->colval;

    // compute idf scores
    idfs = docs->cwgts = da_vmalloc(ncols, "da_csr_Scale: cscale");
    for (i=0; i<ncols; ++i)
        idfs[i] = (colptr[i+1] > colptr[i] ?
                log(1.0*nrows/(colptr[i+1]-colptr[i])) : 0.0);

    // compute row norms
    norms = docs->rnorms = da_vmalloc(nrows, "docs->rnorms");
    for(i=0; i < nrows; ++i){
        for(norm=0.0, j=rowptr[i]; j<rowptr[i+1]; ++j){
            idf = idfs[rowind[j]];
            norm += rowval[j] * rowval[j] * idf * idf;
        }
        norms[i] = norm > 0 ? sqrt(norm) : 0.0;
    }

    // create block index
    ids   = da_umalloc(nnz, "bmmFindNeighbors: ids");
    freqs = da_umalloc(nnz, "bmmFindNeighbors: freqs");
    for(i=0; i < nnz; ++i){
        ids[i] = (uint) colind[i];
        freqs[i] = (uint) colval[i];
    }
    idx = createBlockIndexC(docs, idfs, norms, ids, freqs);

    // now can free docs col structure if necessary
    da_csr_FreeBase(docs, DA_COL);

    // find maximum row/query length
    for(mrlen=0, i=0; i < nrows; ++i)
        if(rowptr[i+1]-rowptr[i] > mrlen)
            mrlen = rowptr[i+1]-rowptr[i];

    da_stopwctimer(params->timer_7); // indexing time

    da_startwctimer(params->timer_5); // memory allocation time
    // allocate space for result
    neighbors = da_csr_Create();
    neighbors->nrows = neighbors->ncols = docs->nrows;
    params->nghnnz = params->k * docs->nrows;
    nptr = neighbors->rowptr = da_pmalloc(docs->nrows + 1, "msFindNeighbors: neighbors->rowptr");
    nind = neighbors->rowind = da_imalloc(params->nghnnz, "msFindNeighbors: neighbors->rowind");
    nval = neighbors->rowval = da_vmalloc(params->nghnnz, "msFindNeighbors: neighbors->rowval");
    nptr[0] = 0;
    params->neighbors = neighbors;
    // allocate search data structures
    qptrs = (da_bil_t**) da_malloc(mrlen * sizeof(da_bil_t*), "qptrs");
    knn = da_iapqCreate(params->k, nrows);
    pmscore = da_vsmalloc(mrlen+1, 0.0, "pmscore");
    bmscore = da_vsmalloc(mrlen, 0.0, "bmscore");
    da_stopwctimer(params->timer_5); // memory allocation time

    // for each x \in V do
    if(params->verbosity > 0)
        printf("Progress Indicator: ");
    fflush(stdout);
    i=0;

    if(params->cr){
        sprintf(params->crfname, "bmmc-%s-%s-%d.cr", da_getDataset(params),
            da_getStringKey(sim_options, params->sim), params->k);
        sprintf(params->crfname2, "%s0", params->crfname);
        if(da_fexists(params->crfname)){
            da_crRead(params, NULL, NULL, &rid);
            i = rid;
            pct = 100*(i+1.0)/(float)nrows;
        }
    }

    for(; i < nrows; i++){
        da_iapqReset(knn);
        da_bmmcGetSimilarRows(i, docs, knn, idx, ids, qptrs, pmscore, bmscore, params);

        // transfer knn to graph
        for(k=nptr[i], j=0; j < knn->nnodes; ++j){
            nind[k]   = knn->heap[j].key;
            nval[k++] = knn->heap[j].val;
        }
        nptr[i+1] = k;
        params->nSimPairs += knn->nnodes;

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
    da_stopwctimer(params->timer_3); // find neighbors time

    /* finalize search */
    if(params->fpout){
        if(params->verbosity > 0){
            printf("\nSorting neighbors in decreasing similarity value order...\n");
            fflush(stdout);
        }
        da_csr_SortValues(neighbors, DA_ROW, DA_SORT_D);

        if(params->verbosity > 0)
            printf("Writing neighborhood matrix to %s.\n", params->oFile);
        da_csr_Write(params->neighbors, params->oFile, params->fmtWrite, 1, 1);
    }

    /* verify results */
    if (params->vFile)
        verify_knng_results(neighbors, params->vFile, params->verbosity);

    freeBlockIndexC(&idx, ncols);
    da_iapqDestroy(knn);
    da_free((void**)&qptrs, &pmscore, &bmscore, &ids, &freqs, LTERM);

}


void da_bmmcGetSimilarRows(uint rid, da_csr_t *docs, da_iapq_t *knn, da_bil_t **idx, uint *idz,
        da_bil_t **qptrs, val_t *pmscore, val_t *bmscore, params_t *params)
{
    ssize_t i, j, rlen, tid, nnon;
    uint ness, nmdid, cid;
    float msc, sim, sim2, t, idf, sc;
    char check_ess;
    float *bsmax, *lsmax;
    val_t rnorm, cnorm;
    ptr_t *rowptr, *lptr, *dptr;
    idx_t *rowind;
    val_t *rowval, *qfreq, *norms;
    uint *bxp, *flag, *ids, *freqs, *bdmax;
    da_bil_t *qptr, **eptrs;

    rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
    rowind  = docs->rowind; // colid for each nnz
    rowval  = docs->rowval; // val for each nnz
    norms   = docs->rnorms; // norms for all vectors
    qfreq   = rowval + rowptr[rid]; // frequencies for query doc
    rnorm   = norms[rid];   // L2 norm of the query vector
    check_ess = 1;

    // Reset the knn heap
    da_iapqReset(knn);

    // gather the columns for the current query
    rlen = rowptr[rid+1]-rowptr[rid];
    if(rlen < 1)
        return;

    for(i=0, j=rowptr[rid]; i < rlen; ++i, ++j){
        qptr = qptrs[i] = idx[rowind[j]];
        qptr->qfr = qfreq[i];
        // reset counters
        for(qptr->i=0; qptr->i < qptr->nbl-1 && *(qptr->dids + qptr->bptr[qptr->i]) == 0; ++qptr->i);  /* first non-empty block */
        qptr->did = decompressBlock(qptr);
    }

    // initial sorting of lists by max score
    da_bilsorti(rlen, qptrs);
    eptrs = qptrs; // essential lists
    ness  = rlen;  // number of essential lists
    nnon  = rlen - ness;

    // compute prefix max scores for query
    pmscore[1] = COMPUTE_SCORE(qptrs[0]->qfr, qptrs[0]->idf, rnorm) * qptrs[0]->maxscore;
    for(i=1; i < rlen; ++i)
        pmscore[i+1] = pmscore[i] + COMPUTE_SCORE(qptrs[i]->qfr, qptrs[i]->idf, rnorm) * qptrs[i]->maxscore;


    // get first k docs and compute their similarities
    while(knn->nnodes < knn->maxsize){
        cid = nextDocId(eptrs, ness);

        if(cid == UINT_MAX)
            return; // out of lists to check
        if(cid == rid){  // no self-similarity
            listsToDocId(eptrs, ness, cid+1);
            continue;
        }
        sim = computeSim(rid, cid, rnorm, norms[cid], ness, eptrs);
        params->nDotProducts++;
        params->nCandidates++;
        da_iapqInsert(knn, cid, sim);
    }

    // process remaining items while shrinking essential columns based on updated threshold
    do {
        // update essential list boundary if necessary
         t = da_iapqSeeTopVal(knn);
         if(check_ess){
             while(ness > 1 && pmscore[rlen-ness+2] < t){
                 eptrs++;
                 ness--;
             }
             nnon = rlen - ness;  // number of non-essential lists
             check_ess = 0;
         }

         // get next doc and compute partial similarity
         cid = nextDocId(eptrs, ness);

         if(cid == UINT_MAX)
             return; // out of lists to check
         if(cid == rid){  // no self-similarity
             listsToDocId(eptrs, ness, cid+1);
             continue;
         }

         cnorm = norms[cid];
         params->nCandidates++;

         // compute the block maxscore based partial similarity of the essential lists
         for(sim=0.0, i=0; i < ness; ++i)
             sim += COMPUTE_SCORE(eptrs[i]->qfr, eptrs[i]->idf, rnorm) * getMaxScore(eptrs[i], cid);

         if(sim + pmscore[nnon] < t){
             listsToDocId(eptrs, ness, nextLiveBlock(eptrs, ness, cid));
             goto nextcid;
         }

         // replace prefix score by block max based prefix score and check bound again
         // also expand set of lists to check for next live block
         for(sim2=0.0, i=nnon-1; i >= 0; --i){
             bmscore[i] = COMPUTE_SCORE(qptrs[i]->qfr, qptrs[i]->idf, rnorm) * getMaxScore(qptrs[i], cid);
             sim2 += bmscore[i];

             if(pmscore[i] + sim + sim2 < t){
                 listsToDocId(eptrs, ness, nextLiveBlock(qptrs+i, rlen-i, cid));
                 goto nextcid;
             }
         }

         // replace suffix score by exact similarity
         sim = computeSim(rid, cid, rnorm, cnorm, ness, eptrs);

         // prune if partial similarity plus prefix sum is less than threshold
         if(nnon > 0 && sim + sim2 < t)
             continue;

         // traverse non-essential lists in reverse order and prune as appropriate
         for(i=nnon-1; i >= 0; --i){
             qptr = qptrs[i];
             if(goToDocId(qptr, cid) == cid){
                 // replace block-based score estimate for this feature with exact one
                 sim2 -= bmscore[i];
                 sim += COMPUTE_SCORE(qptr->qfr, qptr->idf, rnorm) * COMPUTE_SCORE(getFreq(qptr), qptr->idf, cnorm);
             }

             if(i > 0 && sim + sim2 < t){
                 listsToDocId(eptrs, ness, cid+1);
                 goto nextcid;
             }
         }
         params->nDotProducts++;

         if(sim > t){
             da_iapqInsert(knn, cid, sim);
             check_ess = 1;
         }

         nextcid:
             continue;

    } while(1);

}




/**
 * Compress data of size sz, storing <prefix> and flags within the stream.
 * Store:
 *  - length of ids data (including header) - a.k.a. where freq data starts from the <where> pointer
 *  - number of de-compressed items
 *  - prefix of ids data
 *  - ids compressed data
 *     - # COMPR_BS-sized blocks in the list | flag of list 1 (in one uint)
 *     - flag of list 2 | flag of list 3 (in one uint)
 *     ...
 *     - compressed data of the COMPR_BS-sized blocks
 *  - frequencies compressed data
 *     - # COMPR_BS-sized blocks in the list | flag of list 1 (in one uint)
 *     - flag of list 2 | flag of list 3 (in one uint)
 *     ...
 *     - compressed data of the COMPR_BS-sized blocks
 */
uint* compressList(uint *where, uint *ids, uint *freq, const uint sz, const uint prefix)
{
    uint i, j, k, l, nbl, fl2;
    uint *cmp, *flg;
    int fl1;
    uint cdata[COMPR_BS];

    if(sz == 0){
        *(where) = sz;
        return (where + 1);
    }

    nbl = ceil((float)sz/COMPR_BS); /* number of blocks needed to compress using blocks of size COMPR_BS */
    cmp = where + 3 + (nbl+2)/2; /* advance compression pointer to beginning of ids data */
    *(where) = sz;
    *(where+2) = prefix;
    flg = where + 3;

    /* compress and store ids data */
    for(flg[0] = nbl, k=0, i=0; i < nbl; ++i){
        if(i == nbl-1){
            /* need to clear compression block as block to be compressed is not full*/
            memset(cdata, 0, COMPR_BS * sizeof(uint));
            da_ucopy(sz-k, ids+k, cdata);
        } else
            da_ucopy(COMPR_BS, ids+k, cdata);
        /* compress cdata */
        fl1 = -1;
        for (l = 0; fl1 < 0; l++) {
            fl1 = pack_encode(&cmp, cdata, l);
        }
        /* store the flag */
        if(i % 2 == 0){
            fl2 = flg[i/2];
            flg[i/2] =  fl2 | (((unsigned int)(fl1)) << 16);
        } else {
            flg[(i+1)/2] = (unsigned int) fl1;
        }

        k += COMPR_BS;
    }

    /* store length of ids data & compression flags */
    *(where+1) = cmp - where;

    /* keep track of where frequencies data starts and leave room for flags */
    flg = cmp;
    cmp += (nbl+2)/2;

    /* compress and store frequency data */
    for(flg[0] = nbl, k=0, i=0; i < nbl; ++i){
        if(i == nbl-1){
            /* need to clear compression block as block to be compressed is not full*/
            memset(cdata, 0, COMPR_BS * sizeof(uint));
            da_ucopy(sz-k, freq+k, cdata);
        } else
            da_ucopy(COMPR_BS, freq+k, cdata);
        /* compress cdata */
        fl1 = -1;
        for (l = 0; fl1 < 0; l++) {
            fl1 = pack_encode(&cmp, cdata, l);
        }
        /* store the flag */
        if(i % 2 == 0){
            fl2 = flg[i/2];
            flg[i/2] =  fl2 | (((unsigned int)(fl1)) << 16);
        } else {
            flg[(i+1)/2] = (unsigned int) fl1;
        }

        k += COMPR_BS;
    }

    return cmp;
}


/**
 * De-compress ids or frequencies or both from a compressed list
 * Will de-compress if the pointer to either ids or freq is non-null
 */
uint decompressList(uint *what, uint *ids, uint *freq)
{
    uint i, j, k, n, prefix, lenf, nbl, fl, sz;
    uint *dcmp, *start, *flg;
    uint cdata[COMPR_BS];

    sz     = what[0];  /* number of elements in de-compressed list */
    if(sz == 0){
        return sz;
    }
    lenf   = what[1];  /* length of ids compressed data (including header) */
    prefix = what[2];  /* ids prefix */

    if(ids != NULL && sz){
        flg = dcmp = what + 3; /* nbl & flags start here */
        nbl = (*(dcmp)) & 65535;
        dcmp += (nbl+2)/2;
        for(i=0; i < nbl; ++i){
            /* get the flag */
            if(i % 2 == 0){
                fl = flg[i/2] >> 16;
            } else {
                fl = flg[(i+1)/2] & 65535;
            }
            dcmp = pack_decode(ids + i*COMPR_BS, dcmp, fl, 0, 0);
        }
        /**
         * Note: First element is exact. The remaining elements in the list are summands.
         *      In other words, the true 5th element (assuming list size is >= 5) can be
         *      found by performing a prefix sum of the first 5 elements in the de-compressed list.
         */
        ids[0] = prefix;
    }

    if(freq != NULL && sz){
        flg = dcmp = what + lenf; /* nbl & flags start here */
        nbl = (*(dcmp)) & 65535;
        dcmp += (nbl+2)/2;
        for(i=0; i < nbl; ++i){
            /* get the flag */
            if(i % 2 == 0){
                fl = flg[i/2] >> 16;
            } else {
                fl = flg[(i+1)/2] & 65535;
            }
            dcmp = pack_decode(freq + i*COMPR_BS, dcmp, fl, 0, 0);
        }
        /* Note: frequency prefix of 1 must be added on the fly */
    }

    return sz;

}

/**
 * Un-compress the block and return the first doc id in the block
 */
inline uint decompressBlock(da_bil_t *eptr)
{
    if(eptr->i >= eptr->nbl)
        return UINT_MAX;

    eptr->j = 0;
    eptr->cfl = 0;
    eptr->sz = decompressList(eptr->dids + eptr->bptr[eptr->i], eptr->ids, NULL);
    if(eptr->sz){
        return eptr->ids[0];
    }
    return UINT_MAX;
}


/**
 * Un-compress the block and return the first doc id in the block
 */
inline uint decompressBlockNum(da_bil_t *eptr, uint bid)
{
    if(bid >= eptr->nbl)
        return UINT_MAX;

    eptr->j = 0;
    eptr->cfl = 0;
    eptr->sz = decompressList(eptr->dids + eptr->bptr[bid], eptr->ids, NULL);
    if(eptr->sz){
        return eptr->ids[0];
    }
    return UINT_MAX;
}



/**
 * Retrieve current frequency from list
 */
inline uint getFreq(da_bil_t *eptr)
{
    if(!eptr->cfl){
        decompressList(eptr->dids + eptr->bptr[eptr->i], NULL, eptr->freqs);
        eptr->cfl = 1;
    }
    return eptr->freqs[eptr->j] + 1;
}



/**
 * Get max score of a block that did may be in
 */
inline float getMaxScore(da_bil_t *eptr, uint did)
{
    uint bid = DBL(did, eptr->bxp);  /* block id */
    return bid < eptr->nbl ? eptr->bsmax[bid] : 0.0;
}


/**
 * Get max score of a block that did may be in
 */
inline float getMaxDid(da_bil_t *eptr, uint did)
{
    uint bid = DBL(did, eptr->bxp);  /* block id */
    return bid < eptr->nbl ? eptr->bdmax[bid] : UINT_MAX;
}


/**
 * Retrieve current docId from list
 */
inline uint getDocId(da_bil_t *eptr)
{
    return eptr->did;
}


/**
 * Advance list pointers to the block that may contain docId & position on doc if possible
 */
inline uint goToDocId(da_bil_t *eptr, uint did)
{
    if(did <= eptr->did)
        return eptr->did;

    uint i;

    i = DBL(did, eptr->bxp);
    if(i >= eptr->nbl){ /* go to last block, last element */
        eptr->i = eptr->nbl-1;
        eptr->did = UINT_MAX;
        return UINT_MAX;
    }

    if(i != eptr->i){
        eptr->i = i;
        for( ; eptr->i < eptr->nbl-1 && *(eptr->dids + eptr->bptr[eptr->i]) == 0; ++eptr->i);  /* first non-empty block */
        eptr->did = decompressBlock(eptr);
    }

    while(eptr->did < did && eptr->j < eptr->sz){
        eptr->did += eptr->ids[++eptr->j];
    }

    if(eptr->j >= eptr->sz){ /* ran out of items in block */
        eptr->i++;
        for( ; eptr->i < eptr->nbl-1 && *(eptr->dids + eptr->bptr[eptr->i]) == 0; ++eptr->i);  /* first non-empty block */
        eptr->did = decompressBlock(eptr);
    }

    return eptr->did;

}



/**
 * Get the min max doc id of the lists
 */
inline uint nextLiveBlock(da_bil_t **eptrs, idx_t sz, uint did)
{
    uint i, m, k;

    // find max did in lists
    m = NEXTBLD(did, eptrs[0]->bxp) - 1;
    for(k=0; k < sz; ++k){
        i = NEXTBLD(did, eptrs[k]->bxp) - 1;
        if(i < m)
            m = i;
    }

    // get to next block
    if(m < UINT_MAX)
        m++;

    if(m < did)
        da_errexit("m < did, %u < %u", m, did);

    return m;
}


/**
 * Find the smallest doc id in all the lists
 * Note: ids and freqs are global pointers over all data
 */
inline uint nextDocId(da_bil_t **eptrs, idx_t sz)
{
    uint i, did;
    da_bil_t *eptr;

    // find the smallest max did within the current blocks in each of the lists
    for(did=UINT_MAX, i=0; i < sz; ++i){
        eptr = eptrs[i];
        if(eptr->did < did){  // ran out of blocks
            did = eptr->did;
        }
    }

    return did;
}



/**
 * Advance all list pointers to docId
 */
inline void listsToDocId(da_bil_t **eptrs, idx_t sz, uint did)
{
    uint i;
    for(i=0; i < sz; ++i)
        goToDocId(eptrs[i], did);

}


/**
 * Advance list pointer to next doc in block or next block
 */
inline uint nextDoc(da_bil_t *eptr)
{
    if(eptr->did == UINT_MAX) /* end of list reached */
        return UINT_MAX;

    if(eptr->j < eptr->sz-1){
        eptr->j++;
        eptr->did += eptr->ids[eptr->j];
        return eptr->did;
    }

    if(eptr->i < eptr->nbl-1){
        eptr->i++;
        for( ; eptr->i < eptr->nbl-1 && *(eptr->dids + eptr->bptr[eptr->i]) == 0; ++eptr->i);  /* first non-empty block */
        eptr->did = decompressBlock(eptr);
        return eptr->did;
    }

    eptr->did = UINT_MAX;
    return UINT_MAX;
}


inline float computeSim(idx_t rid, idx_t did, val_t rnorm, val_t cnorm, idx_t sz, da_bil_t **eptrs)
{
    idx_t i;
    float sim;
    da_bil_t *eptr;

    for(sim=0.0, i=0; i < sz; ++i){
        eptr = eptrs[i];
        if(goToDocId(eptr, did) == did){
            sim += COMPUTE_SCORE(eptr->qfr, eptr->idf, rnorm) * COMPUTE_SCORE(getFreq(eptr), eptr->idf, cnorm);
            nextDoc(eptr);
        }
    }

    return sim;
}



/**
 * Create a block index from the document term matrix
 */
da_bil_t ** createBlockIndexC(da_csr_t *docs, val_t *idfs, val_t *norms, uint *ids, uint *freqs)
{
    ssize_t i, j, k, l, n, nrows, ncols, nnz;
    uint did, mbdid, df, bx, bs, mbl, mcl, cl;
    float mbs, mls, sc;
    ptr_t *colptr;
    idx_t *colind;
    uint *bxp, *nbl, *cnt, *dli, *dlf, *cid, *cfr, *start;
    da_bil_t** idx, *list;

    nrows  = docs->nrows;
    ncols  = docs->ncols;
    nnz    = docs->rowptr[nrows];

    colptr = docs->colptr; // index in colind/colval/colids where each column starts
    colind = docs->colind; // rowid for each nnz in inv index

    idx = (da_bil_t**) da_malloc(ncols * sizeof(da_bil_t*), "createBlockIndex: idx");
    for(mcl=0, i=0; i < ncols; ++i){
        idx[i] = (da_bil_t*) da_malloc(sizeof(da_bil_t), "createBlockIndex: idx[i");
        memset(idx[i], 0, sizeof(da_bil_t));
        if(colptr[i+1]-colptr[i] > mcl)
            mcl = colptr[i+1]-colptr[i];
    }

    bxp = da_usmalloc(ncols, 0, "createBlockIndex: bxp");
    nbl = da_usmalloc(ncols, 0, "createBlockIndex: nbl");
    dli = da_umalloc(mcl, "createBlockIndex: dli"); // delta list for ids
    dlf = da_umalloc(mcl, "createBlockIndex: dlf"); // delta list for frequencies

    // find the max column length, and the exponent and number of blocks for each list
    for(i=0; i < ncols; ++i){
        if(colptr[i+1] == colptr[i])
            continue;
        did = colind[colptr[i+1]-1]; // max doc id in list
        j = 0;
        while(did > 1<<j)
            j++;
        bxp[i] = BLSIZES[j];
        nbl[i] = DBL(did, bxp[i]) + 1;
    }
    cnt = da_umalloc(da_umax(ncols, nbl), "createBlockIndex: cnt");

    // transfer data to block index
    for(nnz=0, i=0; i < ncols; ++i){
        if(colptr[i+1] == colptr[i])
            continue;
        bx = bxp[i]; // block exponent for this list

        // create dids and freqs deltas
        cl = colptr[i+1] - colptr[i];
        for(j=0, k=colptr[i]; j < cl; ++j, ++k){
            dli[j] = ids[k];
            dlf[j] = freqs[k] - 1;
        }
        for (j = cl-1; j > 0; --j)
            dli[j] -= dli[j-1];

        list      = idx[i];
        list->tid = i;
        list->idf = idfs[i];
        list->bxp = bx;
        list->nbl = n = nbl[i];

        list->bptr   = da_pmalloc(n+1, "createBlockIndex: list->bptr");
        list->bdmax  = da_usmalloc(n, UINT_MAX, "createBlockIndex: list->bdmax");
        list->bsmax  = da_fsmalloc(n, 0.0, "createBlockIndex: list->bsmax");
        list->dids   = da_usmalloc((COMPR_BS+1)*n, 0, "createBlockIndex: list->dids");  // this is where compressed data will be stored
        cfr          = list->dids;

        // count the number of docs in each block in this list
        da_usetzero(n, cnt);
        for(j=colptr[i]; j < colptr[i+1]; ++j)
            cnt[ DBL((uint)colind[j], bx) ]++;

        mbl = da_umax(n, cnt);  /* max items in a block */
        list->ids    = da_usmalloc(mbl+COMPR_BS, 0, "createBlockIndex: list->ids"); // this is for de-compressed data
        list->freqs  = da_usmalloc(mbl+COMPR_BS, 0, "createBlockIndex: list->freqs");

        // compress the block data
        for(list->bptr[0]=0, k=0, j=0; j < n; ++j){
            cfr = compressList(cfr, dli+k, dlf+k, cnt[j], *(ids+nnz+k));
            list->bptr[j+1] = cfr - list->dids;
            k += cnt[j];
        }
        list->dids = da_urealloc(list->dids, list->bptr[n], "Shrinking compressed list");

        /* find max dids and scores for blocks */
        for(j=colptr[i]; j < colptr[i+1]; ++j){
            did = ids[nnz];
            df = freqs[nnz++];
            sc = COMPUTE_SCORE(df, idfs[i], norms[did]);
            k = DBL(did, bx);
            if(did > list->bdmax[k])
                list->bdmax[k] = did;
            if(sc > list->bsmax[k])
                list->bsmax[k] = sc;
        }
        list->maxscore = da_fmax(n, list->bsmax);
        list->did = list->dids[0]; /* starting doc in first block */
        list->max = list->bsmax[0]; /* max of first block*/

        /* fix max dids for 0-length blocks -- needed for de-compression */
        for(j=1; j < n; ++j){
            if(list->bdmax[j] < list->bdmax[j-1])
                list->bdmax[j] = list->bdmax[j-1];
        }
    }

    da_free((void**)&nbl, &bxp, &cnt, &dli, &dlf, LTERM);

    return idx;
}

void freeBlockIndexC(da_bil_t*** idxp, idx_t sz)
{
    size_t i;
    da_bil_t **idx;

    idx = *idxp;

    for(i=0; i < sz; ++i)
        da_free((void**)&idx[i]->bsmax, &idx[i]->bdmax, &idx[i]->bptr,
                &idx[i]->dids, &idx[i]->ids, &idx[i]->freqs,
                &idx[i], LTERM);

    da_free((void**)idxp, LTERM);
}

