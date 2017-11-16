/*!
 \file  maxscore.c
 \brief This file contains Maxscore related functions.
         Maxscore is adapted to work with the cosine function, following the method in the paper below.
         Then, a k-nn query is executed with Maxscore for each document in the dataset to obtain the final
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

// forward declarations
inline idx_t nextDocId(da_il_t **eptrs, idx_t sz);
inline void goToDocId(da_il_t **eptrs, idx_t sz, idx_t did);
inline void goToDocIdList(da_il_t *eptr, idx_t did);
inline float computeSim(idx_t rid, idx_t did, val_t rnorm, val_t cnorm, idx_t sz, da_il_t **eptrs);
void da_msGetSimilarRows(idx_t rid, da_csr_t *docs, da_iapq_t *knn, da_il_t **qptrs,
        val_t *pmscore, params_t *params);


/**
 * Main entry point to KnnIdxJoin.
 */
void msFindNeighbors(params_t *params)
{

    ssize_t i, j, k, nneighbs;
    size_t rid, nnz;
    idx_t nrows, ncols, progressInd, pct, mrlen;
    ptr_t *rowptr, *colptr, *nptr;
    idx_t *rowind, *colind, *nind, *cfreq;
    val_t *rowval, *colval, *nval, *mscore, *idfs, *norms, *pmscore;
    da_csr_t *docs, *neighbors=NULL;
    da_il_t ** qptrs = NULL;
    da_iapq_t *knn = NULL;
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
    ncols   = docs->ncols;  // num cols

    // create inverted index and sort inv lists in doc id order
    da_csr_CreateIndex(docs, DA_COL);
    da_csr_SortIndices(docs, DA_COL);
    colptr  = docs->colptr; // index in colind/colval/colids where each column starts
    colind  = docs->colind; // rowid for each nnz in inv index
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


    // transfer frequencies to colids & later free colval
    cfreq   = docs->colids = da_imalloc(nnz, "docs->colids");
    for(i=0; i < nnz; ++i)
        cfreq[i] = (idx_t) colval[i];

    // find maxscore for each column - cvols for max column scores
    mscore = docs->cvols = da_vsmalloc(ncols, 0.0, "docs->cvols");
    for(i=0; i < ncols; ++i){
        for(msc=0.0, j=colptr[i]; j < colptr[i+1]; ++j){
            sc = (colval[j] * idfs[i]) / norms[colind[j]];
            if (sc > msc)
                msc = sc;
        }
        mscore[i] = msc;
    }
    // can free docs->colval if necessary
    da_free((void**)&docs->colval, LTERM);

    // find maximum row/query length
    for(mrlen=0, i=0; i < nrows; ++i)
        if(rowptr[i+1]-rowptr[i] > mrlen)
            mrlen = rowptr[i+1]-rowptr[i];

    da_stopwctimer(params->timer_7); // indexing time


    // allocate memory
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
    qptrs = (da_il_t**) da_malloc(mrlen * sizeof(da_il_t*), "qptrs");
    for(i=0; i < mrlen; ++i){
        qptrs[i] = (da_il_t*) da_malloc(sizeof(da_il_t), "qptrs[i]");
        memset(qptrs[i], 0, sizeof(da_il_t));
    }
    knn = da_iapqCreate(params->k, nrows);
    pmscore = da_vmalloc(mrlen, "pmscore");
    da_stopwctimer(params->timer_5); // memory allocation time


    // for each x \in V do
    if(params->verbosity > 0)
        printf("Progress Indicator: ");
    i=0;

    if(params->cr){
        sprintf(params->crfname, "msc-%s-%s-%d.cr", da_getDataset(params),
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
        da_msGetSimilarRows(i, docs, knn, qptrs, pmscore, params);

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

    da_iapqDestroy(knn);
    for(i=0; i < mrlen; ++i)
        da_free((void**)&qptrs[i], LTERM);
    da_free((void**)&qptrs, &pmscore, LTERM);

}

void da_msGetSimilarRows(idx_t rid, da_csr_t *docs, da_iapq_t *knn, da_il_t **qptrs,
        val_t *pmscore, params_t *params)
{
    ssize_t i, j, rlen, cid, tid, nnon;
    uint ness, nmdid;
    float msc, sim, t, idf;
    char check_ess;
    val_t rnorm, cnorm;
    ptr_t *rowptr, *colptr;
    idx_t *rowind, *colind, *cfreq;
    val_t *rowval, *qfreq, *mscore, *idfs, *norms;
    da_il_t *qptr, **eptrs;

    rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
    rowind  = docs->rowind; // colid for each nnz
    rowval  = docs->rowval; // val for each nnz
    norms   = docs->rnorms; // norms for all vectors
    idfs    = docs->cwgts;  // idfs for all columns
    qfreq   = rowval + rowptr[rid]; // frequencies for query doc
    rnorm   = norms[rid];   // L2 norm of the query vector
    check_ess = 1;

    colptr  = docs->colptr; // index in colind/colval/colids where each column starts
    colind  = docs->colind; // rowid for each nnz in inv index
    cfreq   = docs->colids;
    mscore  = docs->cvols;

    // Reset the knn heap
    da_iapqReset(knn);

    // gather the columns for the current query
    rlen = rowptr[rid+1]-rowptr[rid];
    if(rlen < 1)
        return;

    for(i=0, j=rowptr[rid]; i < rlen; ++i, ++j){
        qptr = qptrs[i];
        tid = qptr->tid = rowind[j];
        qptr->qfr = qfreq[i];
        qptr->idf = idfs[tid];
        qptr->maxscore = mscore[tid];
        qptr->size = colptr[tid+1] - colptr[tid];
        qptr->freqs  = cfreq + colptr[tid];
        qptr->ids    = colind + colptr[tid];
    }

    // initial sorting of lists by max score
    da_ilsorti(rlen, qptrs);
    eptrs = qptrs; // essential lists
    ness  = rlen;  // number of essential lists
    nnon  = rlen - ness;

    // compute prefix max scores for query
    pmscore[0] = COMPUTE_SCORE(qptrs[0]->qfr, qptrs[0]->idf, rnorm) * qptrs[0]->maxscore;
    for(i=1; i < rlen; ++i)
        pmscore[i] = pmscore[i-1] + COMPUTE_SCORE(qptrs[i]->qfr, qptrs[i]->idf, rnorm) * qptrs[i]->maxscore;

    // get first k docs and compute their similarities
    while(knn->nnodes < knn->maxsize){
        cid = nextDocId(eptrs, ness);
        if(cid == INT_MAX)
            return; // out of lists to check
        if(cid == rid){  // no self-similarity
            goToDocId(eptrs, ness, cid+1);
            continue;
        }
        goToDocId(eptrs, ness, cid);  // advance all pointers to point to cid
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
            while(ness > 1 && pmscore[rlen-ness+1] < t){
                eptrs++;
                ness--;
            }
            nnon = rlen - ness;  // number of non-essential lists
            check_ess = 0;
        }

        // get next doc and compute partial similarity
        cid = nextDocId(eptrs, ness);
        if(cid == INT_MAX)
            return; // out of lists to check
        if(cid == rid){  // no self-similarity
            goToDocId(eptrs, ness, cid+1);
            continue;
        }
        goToDocId(eptrs, ness, cid);  // advance all pointers in essential lists to point to cid
        cnorm = norms[cid];
        sim = computeSim(rid, cid, rnorm, cnorm, ness, eptrs);
        params->nCandidates++;

        // prune if partial similarity plus prefix sum is less than threshold
        if(nnon > 0 && sim + pmscore[nnon-1] < t)
            continue;

        // traverse non-essential lists in reverse order and prune as appropriate
        for(i=nnon-1; i >= 0; --i){
            qptr = qptrs[i];
            goToDocIdList(qptr, cid);
            if(qptr->size > 0 && *qptr->ids == cid){
                sim += COMPUTE_SCORE(qptr->qfr, qptr->idf, rnorm) * COMPUTE_SCORE(*qptr->freqs, qptr->idf, cnorm);
                qptr->freqs++;
                qptr->ids++;
                qptr->size--;
            }
            if(i > 0 && sim + pmscore[i-1] < t)
                goto nextcid;
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
 * Find the smallest doc id in all the lists
 */
inline idx_t nextDocId(da_il_t **eptrs, idx_t sz)
{
    idx_t i, did;
    for(did=INT_MAX, i=0; i < sz; ++i)
        if(eptrs[i]->size > 0 && eptrs[i]->ids[0] < did)
            did = eptrs[i]->ids[0];

    return did;
}



/**
 * Advance all list pointers to docId
 */
inline void goToDocId(da_il_t **eptrs, idx_t sz, idx_t did)
{
    idx_t i;
    da_il_t *eptr;
    for(i=0; i < sz; ++i){
        eptr = eptrs[i];
        while(eptr->size > 0 && *eptr->ids < did){
            eptr->freqs++;
            eptr->ids++;
            eptr->size--;
        }
    }
}

/**
 * Advance a list pointer to docId
 */
inline void goToDocIdList(da_il_t *eptr, idx_t did)
{
    idx_t i;
    while(eptr->size > 0 && *eptr->ids < did){
        eptr->freqs++;
        eptr->ids++;
        eptr->size--;
    }
}

inline float computeSim(idx_t rid, idx_t did, val_t rnorm, val_t cnorm, idx_t sz, da_il_t **eptrs)
{
    idx_t i;
    float sim;
    da_il_t *eptr;

    for(sim=0.0, i=0; i < sz; ++i){
        eptr = eptrs[i];
        if(eptr->size > 0 && *eptr->ids == did){
            sim += COMPUTE_SCORE(eptr->qfr, eptr->idf, rnorm) * COMPUTE_SCORE(*eptr->freqs, eptr->idf, cnorm);
            eptr->freqs++;
            eptr->ids++;
            eptr->size--;
        }
    }

    return sim;
}

