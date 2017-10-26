/*!
 \file  mk.c
 \brief Implementation of bounded Tanimoto search from "Bounds on Lengths of Real Values Vectors Similar with Regard to the Tanimoto Similarity"
 by Marzena Kryszkiewicz, ACIIDS 2013, Part 1
 + Using Non-Zero Dimensions for the Cosine and Tanimoto Similarity Search

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declaration
val_t mkComputeSimilarity(da_csr_t *mat, idx_t i1, idx_t i2, char simtype);

/**
 * Main entry point to mkjoin.
 */
void mkjFindNeighbors(params_t *params)
{

    ssize_t i, j, k, l, c, nsims, np, psz, nneighbs;
    idx_t nrows, ncols, cid, ncands, progressInd;
    val_t v, s, rnrm, rnrm2, cnrm, rnmin, rnmax;
    double ns, a, b, t;
    ptr_t *rowptr, *colptr;
    idx_t *rowind, *colind, *cands = NULL;
    char *mark = NULL;
    val_t *rowval, *colval, *rnorms = NULL, *sims = NULL;
    val_t simT;
    da_ivkv_t *sort = NULL, *hits = NULL;
    da_csr_t *docs, *neighbors = NULL;

    if(params->sim != DA_SIM_TAN && params->sim != DA_SIM_JAC){
        gk_errexit(SIGERR, "The mkJoin method is only implemented for Tanimoto.\n");
    }

    docs    = params->docs;
    nrows   = docs->nrows;  // num rows
    ncols   = docs->ncols;  // num cols
    rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
    rowind  = docs->rowind; // colid for each nnz
    rowval  = docs->rowval; // val for each nnz
    simT    = params->simT - TANEPS;
    t       = ((double)params->simT)/(1.0+((double)params->simT));
    t       = 1.0 - 4.0*t*t;
    psz     = 0;
#ifdef MKPL
    b       = 1.0 + 1.0/((double)params->simT); // 1 + 1/e
    a       = 0.5 * (b+ sqrt(b*b-4.0));         // alpha
#endif

    progressInd = nrows/20;
    ncands  = 0; // number of considered candidates

    // allocate memory
    gk_startwctimer(params->timer_5); // memory allocation time
    sort   = da_ivkvsmalloc(ncols, (da_ivkv_t) {0, 0.0}, "mkFindNeighbors: sort");
    sims   = da_vsmalloc(nrows, 0, "mkFindNeighbors: sims");
    mark   = da_csmalloc(nrows, S_INIT, "mkFindNeighbors: mark");
    cands  = da_ismalloc(nrows, 0, "mkFindNeighbors: cands");
    gk_stopwctimer(params->timer_5); // memory allocation time

    gk_startwctimer(params->timer_7); // indexing time
    if(params->sim == DA_SIM_TAN && params->norm > 0){
        printf("Warning: setting -norm=0 since -sim=tan.\n");
        params->norm = 0;
    }
    simSearchSetup(params);
    da_csr_SortIndices(docs, DA_ROW);
    da_csr_ComputeSquaredNorms(docs, DA_ROW);
    rnorms = docs->rnorms;
    for(i=0; i < nrows; ++i){
        rnorms[i] = sqrt(rnorms[i]);
    }

    da_csr_Normalize(docs, DA_ROW, 2);
    da_csr_CreateIndex(docs, DA_COL);
    colptr  = docs->colptr;
    colind  = docs->colind;
    colval  = docs->colval;
    gk_stopwctimer(params->timer_7); // indexing time

    if(params->nim)
        neighbors = params->neighbors;

    // for each x \in V do
    gk_startwctimer(params->timer_3); // find neighbors time
    if(params->verbosity > 0)
        printf("Progress Indicator (no. of records): ");
    for(i=0; i < nrows; i++){
        // find relevant dimensions (set J), locate candidates and compute their similarities to the query
        for(k=0, j=rowptr[i]; j < rowptr[i+1]; ++j, ++k ){
            sort[k].key = rowind[j];
            sort[k].val = rowval[j];
        }
        da_ivkvsortd(k, sort);
        rnrm = rnorms[i];
        rnrm2 = rnrm * rnrm;
#ifdef MKPL
        rnmin = rnrm * (1.0/a);
        rnmax = rnrm * a;
#endif
        gk_startwctimer(params->timer_1);
        for(ns=0.0, ncands=0, np=0, j=0; j < k; ++j){
            // find candidates and compute part of their dot-products
            c = sort[j].key;
            v = sort[j].val;
            for(l=colptr[c]; l < colptr[c+1]; ++l){
                cid = colind[l];
                if(cid >= i){
                    // only interested in smaller rows
                    break;
                }
                if(mark[cid] == S_INIT){
                    cands[ncands] = cid;
                    ncands++;
                    mark[cid] = S_START;

#ifdef MKPL
                    // prune candidates based on vector length bound
                    cnrm = rnorms[cid];
                    if(cnrm > rnmax || cnrm < rnmin){
                        mark[cid] = S_PRUNED;
                        np++;
                    }
#endif
                    sims[cid] += v * colval[l];
                } else if(mark[cid] == S_START){
                    sims[cid] += v * colval[l];
                }

            }
            if(ns > t){
                j++;
                break;
            }
            ns += v * v;
        }
        psz += j;  // prefix size
        // complete dot-product computation for candidates
        gk_stopwctimer(params->timer_1);
        gk_startwctimer(params->timer_2);
        for( ; j < k; ++j){
            // find candidates and compute their similarity
            c = sort[j].key;
            v = sort[j].val;
            for(l=colptr[c]; l < colptr[c+1]; ++l){
                cid = colind[l];
                if(cid >= i){
                    // only interested in smaller rows
                    break;
                }
                if(mark[cid] != S_INIT){
                    sims[cid] += v * colval[l];
                }
            }
        }
        gk_stopwctimer(params->timer_2);
        params->nPruneTanLen += np;
        params->nCandidates += ncands;
        params->nDotProducts += ncands - np;

        // scale dot-products and find number of neighbors
        for(nsims=0, j=0; j < ncands; ++j){
            if(mark[cands[j]] == S_START){
                cid = cands[j];
                sims[cid] /= ((rnrm2 + rnorms[cid]*rnorms[cid])/(rnrm*rnorms[cid]) - sims[cid]);
                if(sims[cid] >= simT){
                    nsims++;
                }
            }
        }

        //transfer candidates and clear sims
        params->nSimPairs += nsims;
        if(params->nim){
            // set current num neighbors for this row's search
            neighbors->rowptr[i+1] = neighbors->rowptr[i];
            // grow neighborhood matrix if necessary
            if(neighbors->rowptr[i] + nsims > params->nghnnz){
                params->nghnnz += gk_max(nsims, params->nghnnz/2);
                da_csr_Grow(neighbors, params->nghnnz);
                params->nghinccnt++;
            }
        }
        for(j=0; j < ncands; j++){
            if(mark[cands[j]] == S_START && sims[cands[j]] >= params->simT){
                da_storeSim(params, i, cands[j], sims[cands[j]]);
            }
            sims[cands[j]] = 0.0;
            mark[cands[j]] = S_INIT;
        }

        if ( params->verbosity > 0 && i % progressInd == 0 ) {
            printf("%d%%..", (int) ceil(100*i/(float)nrows));
            fflush(stdout);
        }
    }
    if(params->verbosity > 0)
        printf("100%%\n");
    gk_stopwctimer(params->timer_3); // find neighbors time

    printf("Average prefix size: %.4f\n", (psz / (double) nrows));

    simSearchFinalize(params);
    gk_free((void**)&cands, &sort, &sims, &mark, LTERM);
}



/**
 * Main entry point to mkjoin.
 */
void mkjFindNeighbors2(params_t *params)
{

    ssize_t i, j, ii, jj, k, l, c, nsims, np, psz, nneighbs;
    idx_t nrows, ncols, cid, ncands, progressInd;
    val_t v, rnrm, rnrm2, cnrm, rnmin, rnmax;
    double ns, a, s, b, t;
    ptr_t *rowptr, *colptr;
    idx_t *rowind, *colind, *cands = NULL;
    char *mark = NULL;
    val_t *rowval, *colval, *rnorms = NULL, *sims = NULL;
    val_t simT;
    da_ivkv_t *sort = NULL, *hits = NULL;
    da_csr_t *docs, *neighbors = NULL;

    if(params->sim != DA_SIM_TAN && params->sim != DA_SIM_JAC){
        gk_errexit(SIGERR, "The mkJoin method is only implemented for Tanimoto.\n");
    }

    docs    = params->docs;
    nrows   = docs->nrows;  // num rows
    ncols   = docs->ncols;  // num cols
    rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
    rowind  = docs->rowind; // colid for each nnz
    rowval  = docs->rowval; // val for each nnz
    simT    = params->simT - TANEPS;
    t       = ((double)params->simT)/(1.0+((double)params->simT));
    t       = 1.0 - 4.0*t*t;
    psz     = 0;


    progressInd = nrows/20;
    ncands  = 0; // number of considered candidates

    // allocate memory
    gk_startwctimer(params->timer_5); // memory allocation time
    sort   = da_ivkvsmalloc(ncols, (da_ivkv_t) {0, 0.0}, "mkFindNeighbors: sort");
    sims   = da_vsmalloc(nrows, 0, "mkFindNeighbors: sims");
    mark   = da_csmalloc(nrows, S_INIT, "mkFindNeighbors: mark");
    cands  = da_ismalloc(nrows, 0, "mkFindNeighbors: cands");
    gk_stopwctimer(params->timer_5); // memory allocation time

    gk_startwctimer(params->timer_7); // indexing time
    if(params->sim == DA_SIM_TAN && params->norm > 0){
        printf("Warning: setting -norm=0 since -sim=tan.\n");
        params->norm = 0;
    }
    simSearchSetup(params);
    da_csr_SortIndices(docs, DA_ROW);
    da_csr_ComputeSquaredNorms(docs, DA_ROW);
    rnorms = docs->rnorms;
    for(i=0; i < nrows; ++i){
        rnorms[i] = sqrt(rnorms[i]);
    }

    da_csr_Normalize(docs, DA_ROW, 2);
    da_csr_CreateIndex(docs, DA_COL);
    colptr  = docs->colptr;
    colind  = docs->colind;
    colval  = docs->colval;
    gk_stopwctimer(params->timer_7); // indexing time

    if(params->nim)
        neighbors = params->neighbors;

    // for each x \in V do
    gk_startwctimer(params->timer_3); // find neighbors time
    if(params->verbosity > 0)
        printf("Progress Indicator (no. of records): ");
    for(i=0; i < nrows; i++){
        // find relevant dimensions (set J), locate candidates and compute their similarities to the query
        for(k=0, j=rowptr[i]; j < rowptr[i+1]; ++j, ++k ){
            sort[k].key = rowind[j];
            sort[k].val = rowval[j];
        }
        da_ivkvsortd(k, sort);
        rnrm = rnorms[i];
        rnrm2 = rnrm * rnrm;
        gk_startwctimer(params->timer_1);
        for(ns=0.0, ncands=0, np=0, j=0; j < k; ++j){
            // find candidates and compute part of their dot-products
            c = sort[j].key;
            v = sort[j].val;

            for(l=colptr[c]; l < colptr[c+1]; ++l){
                cid = colind[l];
                if(cid >= i){
                    // only interested in smaller rows
                    break;
                }
                if(mark[cid] == S_INIT){
                    cands[ncands] = cid;
                    ncands++;
                    mark[cid] = S_START;

                    // find values in candidate that do not exist in query
                    s = 0.0;
                    ii = rowptr[i];
                    jj = rowptr[cid];
                    for( ; ii < rowptr[i+1] && jj < rowptr[cid+1]; ){
                        if(rowind[jj] == rowind[ii]){
                            ii++;
                            jj++;
                        } else
                        if(rowind[ii] < rowind[jj]){
                            s += rowval[ii] * rowval[ii];
                            ii++;
                        } else {
                            jj++;
                        }
                    }
                    for( ; ii < rowptr[i+1]; ++ii){
                        s += rowval[ii] * rowval[ii];
                    }
                    b = (1.0 + 1.0/((double)params->simT)) * sqrt(1.0-s); // (1 + 1/e)*sqrt(1-\sum_{i\in L}{u_i^2}), L is the set of non-zero query features which are zero candidate features
                    a = 0.5 * (b + sqrt(b*b-4.0));                       // beta
                    rnmin = rnrm * (1.0/a);
                    rnmax = rnrm * a;

                    // prune candidates based on vector length bound
                    cnrm = rnorms[cid];

                    if(b < 2 || cnrm > rnmax || cnrm < rnmin){
                        mark[cid] = S_PRUNED;
                        np++;
                        continue;
                    }

                    sims[cid] += v * colval[l];
                } else if(mark[cid] == S_START){
                    sims[cid] += v * colval[l];
                }

            }
            if(ns > t){
                j++;
                break;
            }
            ns += v * v;
        }
        psz += j;  // prefix size
        // complete dot-product computation for candidates
        gk_stopwctimer(params->timer_1);
        gk_startwctimer(params->timer_2);
        for( ; j < k; ++j){
            // find candidates and compute their similarity
            c = sort[j].key;
            v = sort[j].val;
            for(l=colptr[c]; l < colptr[c+1]; ++l){
                cid = colind[l];
                if(cid >= i){
                    // only interested in smaller rows
                    break;
                }
                if(mark[cid] != S_INIT){
                    sims[cid] += v * colval[l];
                }
            }
        }
        gk_stopwctimer(params->timer_2);
        params->nPruneTanLen += np;
        params->nCandidates += ncands;
        params->nDotProducts += ncands - np;

        // scale dot-products and find number of neighbors
        for(nsims=0, j=0; j < ncands; ++j){
            if(mark[cands[j]] == S_START){
                cid = cands[j];
                sims[cid] /= ((rnrm2 + rnorms[cid]*rnorms[cid])/(rnrm*rnorms[cid]) - sims[cid]);
                if(sims[cid] >= simT){
                    nsims++;
                }
            }
        }

        //transfer candidates and clear sims
        params->nSimPairs += nsims;
        if(params->nim){
            // set current num neighbors for this row's search
            neighbors->rowptr[i+1] = neighbors->rowptr[i];
            // grow neighborhood matrix if necessary
            if(neighbors->rowptr[i] + nsims > params->nghnnz){
                params->nghnnz += gk_max(nsims, params->nghnnz/2);
                da_csr_Grow(neighbors, params->nghnnz);
                params->nghinccnt++;
            }
        }
        for(j=0; j < ncands; j++){
            if(mark[cands[j]] == S_START && sims[cands[j]] >= params->simT){
                da_storeSim(params, i, cands[j], sims[cands[j]]);
            }
            sims[cands[j]] = 0.0;
            mark[cands[j]] = S_INIT;
        }

        if ( params->verbosity > 0 && i % progressInd == 0 ) {
            printf("%d%%..", (int) ceil(100*i/(float)nrows));
            fflush(stdout);
        }
    }
    if(params->verbosity > 0)
        printf("100%%\n");
    gk_stopwctimer(params->timer_3); // find neighbors time

    printf("Average prefix size: %.4f\n", (psz / (double) nrows));

    simSearchFinalize(params);
    gk_free((void**)&cands, &sort, &sims, &mark, LTERM);
}


val_t mkComputeSimilarity(da_csr_t *mat, idx_t i1, idx_t i2, char simtype)
{
    idx_t nind1, nind2;
    idx_t *ind1, *ind2;
    val_t *val1, *val2, *rnorm, sim, n1, n2;

    if (!mat->rowptr)
        gk_errexit(SIGERR, "Row-based view of the matrix does not exists.\n");
    if (!mat->rnorms)
        gk_errexit(SIGERR, "Row-based squared norms not present.\n");
    nind1 = mat->rowptr[i1+1]-mat->rowptr[i1];
    nind2 = mat->rowptr[i2+1]-mat->rowptr[i2];
    ind1  = mat->rowind + mat->rowptr[i1];
    ind2  = mat->rowind + mat->rowptr[i2];
    val1  = mat->rowval + mat->rowptr[i1];
    val2  = mat->rowval + mat->rowptr[i2];
    n1    = mat->rnorms[i1];
    n2    = mat->rnorms[i2];


    sim = 0.0;
    i1 = i2 = 0;
    while (i1<nind1 && i2<nind2) {
        if (i1 == nind1) {
            i2++;
        }
        else if (i2 == nind2) {
            i1++;
        }
        else if (ind1[i1] < ind2[i2]) {
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
    if (simtype == DA_SIM_COS)
        sim = (n1*n2 > 0.0 && n1*n2 != 1.0 ? sim/sqrt(n1*n2) : 0.0);
    else
        sim = (n1+n2-sim > 0.0 ? sim/(n1+n2-sim) : 0.0);

    return sim;

}
