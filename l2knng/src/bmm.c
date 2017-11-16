/*!
 \file  bmm.c
 \brief This file contains Block-Max Maxscore (BMM) related functions.
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

#define dl_bitsize(a)\
    (8 * sizeof(a))

#if defined(__GNUC__) || defined(__INTEL_COMPILER)
   #define dl_clz(a) \
     ((size_t)( (a) == 0 ? dl_bitsize(a) : \
       (sizeof(a) == sizeof(unsigned long long) ? \
         (size_t)__builtin_clzll((unsigned long long)a) \
         : (size_t)(__builtin_clz((unsigned int)a) - \
           (dl_bitsize(unsigned int) - dl_bitsize(a))) \
       ) \
     ))
 #else
   / cause trouble */
 #endif


// forward declarations
inline float getMaxScore(da_bil_t *eptr, uint did);
inline float getMaxDid(da_bil_t *eptr, uint did);
inline uint getDocId(da_bil_t *eptr);
inline uint getFreq(da_bil_t *eptr);
inline uint goToDocId(da_bil_t *eptr, uint did);
inline void movelistsToDocIdBlock(da_bil_t **eptrs, idx_t sz, uint did);
//inline uint nextLiveBlock(da_bil_t **eptrs, idx_t sz, uint did);
inline uint nextLiveBlock(da_bil_t **eptrs, idx_t sz, idx_t from, idx_t ness, uint did, float t, val_t *norms);
inline uint nextDocId(da_bil_t **eptrs, idx_t sz);
inline void listsToDocId(da_bil_t **eptrs, idx_t sz, uint did);
inline uint nextDoc(da_bil_t *eptr);

inline float computeSim(idx_t rid, idx_t did, val_t rnorm, val_t cnorm, idx_t sz, da_bil_t **eptrs);
void da_bmmGetSimilarRows(uint rid, da_csr_t *docs, da_iapq_t *knn, da_bil_t **idx, uint *idz,
        da_bil_t **qptrs, val_t *pmscore, val_t *bmscore, params_t *params);
da_bil_t** createBlockIndex(da_csr_t *docs, val_t *idfs, val_t *norms, uint *ids, uint *freqs);
void freeBlockIndex(da_bil_t ***idxp, idx_t sz);

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



/**
 * Main entry point to BMM.
 */
void bmmFindNeighbors(params_t *params)
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
    idx = createBlockIndex(docs, idfs, norms, ids, freqs);

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
        sprintf(params->crfname, "bmm-%s-%s-%d.cr", da_getDataset(params),
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
        da_bmmGetSimilarRows(i, docs, knn, idx, ids, qptrs, pmscore, bmscore, params);

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

    freeBlockIndex(&idx, ncols);
    da_iapqDestroy(knn);
    da_free((void**)&qptrs, &pmscore, &bmscore, &ids, &freqs, LTERM);

}


void da_bmmGetSimilarRows(uint rid, da_csr_t *docs, da_iapq_t *knn, da_bil_t **idx, uint *idz,
        da_bil_t **qptrs, val_t *pmscore, val_t *bmscore, params_t *params)
{
    ssize_t i, j, rlen, tid, nnon;
    uint ness, nmdid, cid;
    float msc, sim, sim2, sim3, et, t, idf, sc;
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
        qptr->j = 0;
        qptr->did = qptr->dids[0];
        qptr->i = DBL(qptr->did, qptr->bxp);
        qptr->max = qptr->bdmax[qptr->i];
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
            listsToDocId(eptrs, ness, nextLiveBlock(qptrs, rlen, nnon, ness, cid, t, norms));
            goto nextcid;
        }

        // replace prefix score by block max based prefix score and check bound again
        // also expand set of lists to check for next live block
        for(sim2=0.0, i=nnon-1; i >= 0; --i){
            bmscore[i] = COMPUTE_SCORE(qptrs[i]->qfr, qptrs[i]->idf, rnorm) * getMaxScore(qptrs[i], cid);
            sim2 += bmscore[i];

            if(pmscore[i] + sim + sim2 < t){
                listsToDocId(eptrs, ness, nextLiveBlock(qptrs, rlen, i, ness, cid, t, norms));
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
 * Retrieve current frequency from list
 */
inline uint getFreq(da_bil_t *eptr)
{
    return eptr->dfreqs[eptr->j];
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
        eptr->max = 0.0;
        return UINT_MAX;
    }

    if(i != eptr->i){
        eptr->i = i;
        eptr->j = eptr->bptr[i];
        eptr->did = eptr->dids[eptr->j];
        eptr->max = eptr->bsmax[eptr->i];
    }

    while(eptr->did < did && eptr->j < eptr->bptr[eptr->i+1]){
        eptr->did = eptr->dids[++eptr->j];
    }

    if(eptr->j == eptr->bptr[eptr->i+1]){ /* ran out of items in block */
        eptr->did = eptr->dids[--eptr->j];
        nextDoc(eptr);
    }

    return eptr->did;

}


/**
 * Advance all list pointers to the block that may contain docId
 */
inline void movelistsToDocIdBlock(da_bil_t **eptrs, idx_t sz, uint did)
{
    uint i, k;
    da_bil_t *eptr;
    for(k=0; k < sz; ++k){
        eptr = eptrs[k];
        if(eptr->did == UINT_MAX)  /* have already gotten to end of this list */
            continue;
        i = DBL(did, eptr->bxp);
        if(i >= eptr->nbl){ /* go to last block, last element */
            eptr->i = eptr->nbl-1;
            eptr->did = UINT_MAX;
            eptr->max = 0.0;
            continue;
        }

        if(i != eptr->i){
            eptr->i = i;
            eptr->max = eptr->bsmax[eptr->i];
        }
    }
}



/**
 * Get the min max doc id of the lists
 */
inline uint nextLiveBlock(da_bil_t **eptrs, idx_t sz, idx_t from, idx_t ness, uint did, float t, val_t *norms)
{
    uint i, m, k;

    // find max did in lists
    m = NEXTBLD(did, eptrs[from]->bxp) - 1;
    for(k=from+1; k < sz; ++k){
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

    ptr_t *ptr = eptr->bptr;
    uint i = eptr->i;

    eptr->j++;
    if(eptr->j >= ptr[eptr->nbl]){ /* ran out of items in any of the blocks */
        eptr->did = UINT_MAX;
        eptr->i = eptr->nbl-1;
        eptr->max = 0.0;
        return UINT_MAX;
    }

    eptr->did = eptr->dids[eptr->j];
    eptr->i = DBL(eptr->did, eptr->bxp);
    eptr->max = eptr->bsmax[eptr->i];

    return eptr->did;
}


inline float computeSim(idx_t rid, idx_t did, val_t rnorm, val_t cnorm, idx_t sz, da_bil_t **eptrs)
{
    idx_t i;
    float sim;
    da_bil_t *eptr;

    for(sim=0.0, i=0; i < sz; ++i){
        eptr = eptrs[i];
        if(goToDocId(eptr, did) == did){
            sim += COMPUTE_SCORE(eptr->qfr, eptr->idf, rnorm) * COMPUTE_SCORE(eptr->dfreqs[eptr->j], eptr->idf, cnorm);
            nextDoc(eptr);
        }
    }

    return sim;
}


/**
 * Create a block index from the document term matrix
 */
da_bil_t ** createBlockIndex(da_csr_t *docs, val_t *idfs, val_t *norms, uint *ids, uint *freqs)
{
    ssize_t i, j, k, l, n, mclen, nrows, ncols, nnz;
    uint did, mbdid, df, bx, bs, mbl;
    float mbs, mls, sc;
    ptr_t *colptr;
    idx_t *colind;
    uint *bxp, *nbl, *cnt;
    da_bil_t** idx, *list;

    nrows  = docs->nrows;
    ncols  = docs->ncols;
    nnz    = docs->rowptr[nrows];

    colptr = docs->colptr; // index in colind/colval/colids where each column starts
    colind = docs->colind; // rowid for each nnz in inv index

    idx = (da_bil_t**) da_malloc(ncols * sizeof(da_bil_t*), "createBlockIndex: idx");
    for(i=0; i < ncols; ++i){
        idx[i] = (da_bil_t*) da_malloc(sizeof(da_bil_t), "createBlockIndex: idx[i");
        memset(idx[i], 0, sizeof(da_bil_t));
    }

    bxp = da_usmalloc(ncols, 0, "createBlockIndex: bxp");
    nbl = da_usmalloc(ncols, 0, "createBlockIndex: nbl");

    // find the max column length, and the exponent and number of blocks for each list
    for(mclen=0, i=0; i < ncols; ++i){
        if(colptr[i+1] == colptr[i])
            continue;
        if(colptr[i+1]-colptr[i] > mclen)
            mclen = colptr[i+1]-colptr[i];
        did = colind[colptr[i+1]-1]; // max doc id in list
        j = 0;
        while(did > 1<<j)
            j++;
        bxp[i] = BLSIZES[j];
        nbl[i] = DBL(did, bxp[i]) + 1;
    }
    mbl = da_umax(ncols, nbl);  /* max number of blocks */
    cnt = da_usmalloc(mbl, 0, "createBlockIndex: cnt");

    // transfer data to block index
    for(nnz=0, i=0; i < ncols; ++i){
        if(colptr[i+1] == colptr[i])
            continue;
        bx = bxp[i]; // block exponent for this list

        list = idx[i];
        list->tid = i;
        list->idf = idfs[i];
        list->bxp = bx;
        list->nbl = n = nbl[i];

        list->bptr   = da_pmalloc(n+1, "createBlockIndex: list->bptr");
        list->bdmax  = da_usmalloc(n, UINT_MAX, "createBlockIndex: list->bdmax");
        list->bsmax  = da_fsmalloc(n, 0.0, "createBlockIndex: list->bsz");
        list->dids   = ids + nnz;  // this is where this list starts in the inv index
        list->dfreqs = freqs + nnz;

        // count the number of docs in each block in this list
        da_usetzero(n, cnt);
        for(j=colptr[i]; j < colptr[i+1]; ++j)
            cnt[ DBL((uint)colind[j], bx) ]++;

        // set up the block pointers
        for(list->bptr[0]=0, j=0; j < n; ++j)
            list->bptr[j+1] = list->bptr[j] + cnt[j];

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
    }

    da_free((void**)&nbl, &bxp, &cnt, LTERM);

    return idx;
}

void freeBlockIndex(da_bil_t*** idxp, idx_t sz)
{
    size_t i;
    da_bil_t **idx;

    idx = *idxp;

    for(i=0; i < sz; ++i)
        da_free((void**)&idx[i]->bsmax, &idx[i]->bdmax, &idx[i]->bptr, &idx[i], LTERM);

    da_free((void**)idxp, LTERM);
}

