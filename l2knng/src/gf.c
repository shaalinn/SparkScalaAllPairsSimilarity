/*!
 \file  gf.c
 \brief This file contains GreedyFiltering related functions.
         GreedyFiltering is re-implemented in C based on the Java based implementation provided by the
         algorithm's authors and is described in the following paper.

@incollection{ park2014,
    year = {2014},
    booktitle = {Database Systems for Advanced Applications},
    volume = {8421},
    series = {Lecture Notes in Computer Science},
    title = {Greedy Filtering: A Scalable Algorithm for K-Nearest Neighbor Graph Construction},
    publisher = {Springer-Verlag},
    author = {Park, Youngki and Park, Sungchan and Lee, Sang-goo and Jung, Woosung},
    pages = {327-341},
    language = {English}
}

 \author David C. Anastasiu
 */

#include "includes.h"

/** Forward declaration */
void addNeighbor(idx_t ii, idx_t jj, val_t v, da_csr_t *neigh, idx_t *cnt, idx_t k)
{
    ptr_t *ptr = neigh->rowptr;
    idx_t *ind = neigh->rowind;
    val_t *val = neigh->rowval;

    size_t i, j;
    val_t mv;
    /* less than k neighbors */
    if(cnt[ii] < k){
        for(j=ptr[ii]; j < ptr[ii]+cnt[ii]; ++j)
            if(ind[j] == jj)
                return;
        ind[ptr[ii]+cnt[ii]] = jj;
        val[ptr[ii]+cnt[ii]] = v;
        cnt[ii]++;
        return;
    }
    /* find neighbor with smallest similarity and replace it */
    for(mv = FLT_MAX, i=j=ptr[ii]; j < ptr[ii]+k; ++j){
        if(ind[j] == jj)
            return;
        if(val[j] < mv){
            i = j;
            mv = val[j];
        }
    }
    if(v > mv){
        ind[i] = jj;
        val[i] = v;
    }
}


/**
 * Compute cosine similarity between two rows in the matrix.
 * Matrix is assumed to have normalized rows.
 */
inline float computeSim(da_csr_t *mat, idx_t ii, idx_t jj)
{
    idx_t nind1, nind2, i1, i2;
    idx_t *ind1, *ind2;
    val_t *val1, *val2, stat1, stat2, sim;

    nind1 = mat->rowptr[ii+1]-mat->rowptr[ii];
    nind2 = mat->rowptr[jj+1]-mat->rowptr[jj];
    ind1  = mat->rowind + mat->rowptr[ii];
    ind2  = mat->rowind + mat->rowptr[jj];
    val1  = mat->rowval + mat->rowptr[ii];
    val2  = mat->rowval + mat->rowptr[jj];

    sim = 0.0;
    i1 = i2 = 0;
    while (i1<nind1 && i2<nind2) {
        if (ind1[i1] < ind2[i2]) {
            i1++;
        }
        else if (ind1[i1] > ind2[i2]) {
            i2++;
        }
        else {
            sim   += val1[i1]*val2[i2];
            i1++;
            i2++;
        }
    }

    return sim;

}


/**
 * Main entry point to KnnIdxJoin.
 */
void gfFindNeighbors(params_t *params)
{
    size_t i, ii, j, jj, k, l, nrows, ncols, nnz, mu, p, n, nn, cc;
    idx_t *ind, *cind, *nind, *pcnt, *fcnt, *todo;
    ptr_t *ptr, *cptr, *nptr;
    val_t *nval, v;
    da_csr_t *docs=NULL, *neigh=NULL;

    /* pre-process data */
    preProcessData(params);

    mu     = da_max(params->k, params->mu);
    docs   = params->docs;
    nrows  = docs->nrows;   // num rows
    ncols  = docs->ncols;   // num cols
    nnz    = docs->rowptr[nrows];
    ind    = docs->rowind;
    ptr    = docs->rowptr;

    /** allocate memory */
    todo = da_ismalloc(nrows, 0, "gfFindNeighbors: todo");
    pcnt = da_ismalloc(nrows, 1, "gfFindNeighbors: pcnt");
    fcnt = da_ismalloc(ncols, 0, "gfFindNeighbors: fcnt");
    neigh = da_csr_Create();  /** neighborhood graph */
    neigh->nrows = neigh->ncols = nrows;
    neigh->rowptr = nptr = da_pmalloc(nrows+1, "gfFindNeighbors: nptr");
    neigh->rowind = nind = da_imalloc(params->k * nrows, "gfFindNeighbors: nind");
    neigh->rowval = nval = da_vmalloc(params->k * nrows, "gfFindNeighbors: nind");
    for(nptr[0]=0, i=0; i < nrows; ++i)
        nptr[i+1] = nptr[i] + params->k;
    if(params->verbosity > 0)
        printf("Docs matrix: " PRNT_IDXTYPE " rows, " PRNT_IDXTYPE " cols, "
            PRNT_PTRTYPE " nnz\n", docs->nrows, docs->ncols, docs->rowptr[docs->nrows]);

    /** Prefix selection */
    da_startwctimer(params->timer_3); // overall knn graph construction time

    da_csr_SortValues(docs, DA_ROW, DA_SORT_D);

    p  = 1; /** current prefix */
    nn = nrows;  /** remaining rows to find prefix for */
    for(i=0; i < nrows; ++i){
        todo[i] = i;  /** contains vectors that do not meet the stop condition */
        j = ptr[i];   /** For each vector, set the prefix of the vector to 1 */
        if(j < ptr[i+1])
            fcnt[ind[j]]++;
    }
    do {
        /** Determine which prefixes will be extended */
        for(n=0, i=0; i < nn; ++i){
            k = todo[i];
            for(cc=0, j=ptr[k]; j < ptr[k+1] && j < ptr[k]+p; ++j){
                cc += fcnt[ind[j]]-1;
            }
            if(cc < mu && p < ptr[k+1]-ptr[k]){
                todo[n++] = k;
            }
        }
        nn = n;

        /** Put the vectors into a bucket */
        for(i=0; i < nn; ++i){
            k = todo[i];
            j=ptr[k]+p;
            if(j < ptr[k+1]){
                fcnt[ind[j]]++;
                pcnt[k]++;
            }
        }

        /** Extend the prefixes */
        p++;

    } while(nn > 0);

    /** Create partial inverted index from prefixes */
    params->nCandidates = nnz = da_isum(ncols, fcnt, 1);
    cptr = da_pmalloc(ncols+1, "gfFindNeighbors: cptr");
    cind = da_imalloc(nnz, "gfFindNeighbors: cind");
    for(cptr[0]=0, j=0; j < ncols; ++j)
        cptr[j+1] = cptr[j] + fcnt[j];
    for(i=0; i < nrows; ++i){
        for(j=ptr[i]; j < ptr[i] + pcnt[i] && j < ptr[i+1]; ++j){
            l = ind[j];
            cind[cptr[l]++] = i;
        }
    }
    CSRSHIFT(i, ncols, cptr);

    /** Before calculating similarities, sort each element into ascending order of its dimension number */
    da_csr_SortIndices(docs, DA_ROW);

    /** Calculate the similarity for each candidate pair */
    da_isetzero(nrows, pcnt);  /* reusing pcnt as neighbor counter */
    for(j=0; j < ncols; ++j){
        for(k=cptr[j]; k < cptr[j+1]; ++k){
            ii = cind[k];
            for(l=k+1; l < cptr[j+1]; ++l){
                jj = cind[l];
                v = computeSim(params->docs, ii, jj);
                addNeighbor(ii, jj, v, neigh, pcnt, params->k);
                addNeighbor(jj, ii, v, neigh, pcnt, params->k);
                params->nDotProducts++;
            }
        }
    }
    da_stopwctimer(params->timer_3); // find neighbors time

    /** compress space for rows with less than k neighbors */
    params->nSimPairs = da_isum(nrows, pcnt, 1);
    if(params->nSimPairs < params->k * nrows){
        for(nnz=0, i=0; i < nrows; ++i){
            for(j=nptr[i]; j < nptr[i]+pcnt[i]; ++j){
                nind[nnz]   = nind[j];
                nval[nnz++] = nval[j];
            }
            nptr[i] = nnz;
        }
        CSRSHIFT(i, nrows, nptr);
    }

    /* finalize search */
    if(params->fpout){
        if(params->verbosity > 0){
            printf("\nSorting neighbors in decreasing similarity value order...\n");
            fflush(stdout);
        }
        da_csr_SortValues(neigh, DA_ROW, DA_SORT_D);

        if(params->verbosity > 0)
            printf("Writing neighborhood matrix to %s.\n", params->oFile);
        da_csr_Write(neigh, params->oFile, params->fmtWrite, 1, 1);
    }

    /* verify results */
    if (params->vFile)
        verify_knng_results(neigh, params->vFile, 0);

    /** free memory */
    da_csr_Free(&neigh);
    da_free((void**)&todo, &pcnt, &fcnt, &cptr, &cind, LTERM);

}

