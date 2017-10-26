/*!
 \file  mmjoin.c
 \brief This file contains MMJoin related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
idx_t mmjProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, val_t *lengths, ptr_t *endptr,
		idx_t ncands, idx_t *cands, accum_t *accum, val_t *marker, val_t *marker2, val_t *rwgts, val_t *rfiwgts);
void mmjIndexRow(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdxJ_t *invIdx,
		val_t *rwgts, val_t *rfiwgts, val_t *cwgts);
void mmjFindMatches(idx_t rid, params_t *params, idx_t *cands,
		da_invIdxJ_t *idx, val_t *lengths, da_csr_t *docs, idx_t *rsizes, val_t *rwgts, val_t *rfiwgts, val_t *cwgts,
		ptr_t *endptr, val_t *marker, val_t *marker2);
void mmjReorderDocs(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
		val_t *cwgts, idx_t *rperm, idx_t *cperm);

/**
 * Main entry point to MMJoin.
 */
void mmjFindNeighbors(params_t *params){
	ssize_t i, j, k, nnz;
	idx_t nrows, ncols, mrlen, mclen, progressInd;
	ptr_t *endptr;
	idx_t *cands, *rsizes, *csizes, *rperm, *cperm;
	val_t maxrval, myval;
	val_t *rwgts = NULL, *rfiwgts = NULL, *cwgts = NULL, *marker, *marker2, *lengths = NULL;
	da_invIdxJ_t *invIdx = NULL; // inverted index
	da_csr_t *docs = params->docs;

    if(params->sim == DA_SIM_TAN && params->norm > 0){
        printf("Warning: setting -norm=0 since -sim=tan.\n");
        params->norm = 0;
    }

	simSearchSetup(params);
	nrows  = docs->nrows;   // num rows
	ncols  = docs->ncols;   // num cols
	progressInd = ceil(nrows/(float)20);

	rsizes  = da_imalloc(nrows, "mmjFindNeighbors: rsizes");
	csizes  = da_ismalloc(ncols, 0.0, "mmjFindNeighbors: csizes");
	rperm   = da_imalloc(nrows, "mmjFindNeighbors: rperm");
	cperm   = da_imalloc(ncols, "mmjFindNeighbors: cperm");
	rwgts   = da_vsmalloc(nrows, 0.0, "mmjFindNeighbors: rwgts");
	rfiwgts = da_vmalloc(nrows, "mmjFindNeighbors: rfiwgts");
	cwgts   = da_vsmalloc(ncols, 0.0, "mmjFindNeighbors: cwgts");
	marker  = da_vsmalloc(ncols, 0.0, "mmjFindNeighbors: marker");
	marker2 = da_vsmalloc(ncols, 0.0, "mmjFindNeighbors: marker2");

	// reorder matrix
	mmjReorderDocs(params, &docs, rsizes, csizes, rwgts, cwgts, rperm, cperm);

	nnz     = docs->rowptr[nrows]; // nnz in docs

	// allocate memory for candidates and the inverted indexes -- needs nnz & col counts
	cands   = da_imalloc(nrows, "mmjFindNeighbors: candidates");
	lengths = da_vsmalloc(nnz, 0.0, "mmjFindNeighbors: lengths");
	endptr  = da_pmalloc(nrows, "mmjFindNeighbors: endptr");
	invIdx  = gk_malloc(sizeof(da_invIdxJ_t), "mmjFindNeighbors: invIdx");
	invIdx->ids    = da_ismalloc(nnz, 0, "mmjFindNeighbors: invIdx->ids");
	invIdx->vals   = da_vsmalloc(nnz, 0.0, "mmjFindNeighbors: invIdx->ids");
	invIdx->lens   = da_vsmalloc(nnz, 0.0, "mmjFindNeighbors: invIdx->ids");
	invIdx->starts = da_pmalloc(ncols, "mmjFindNeighbors: invIdx->starts");
	invIdx->ends   = da_pmalloc(ncols, "mmjFindNeighbors: invIdx->ends");
	for(j=0, k=0; j < ncols; j++){
		invIdx->starts[j] = invIdx->ends[j] = k;
		k += csizes[j];
	}
	invIdx->accum  = gk_malloc(nrows*sizeof(accum_t), "mmjFindNeighbors: invIdx->accum");
	for(i=0; i < nrows; i++)
		invIdx->accum[i] = -1;

#ifdef EXTRACOUNTS
	// find the average length of rows & columns
	mrlen = ceil(da_isum(nrows, rsizes, 1) / (val_t) nrows);
	mclen = ceil(da_isum(ncols, csizes, 1) / (val_t) ncols);
	if(params->verbosity > 0)
		printf("mrlen: " PRNT_IDXTYPE ", mclen: " PRNT_IDXTYPE ".\n", mrlen, mclen);
#endif

	// for each x \in V do
	gk_startwctimer(params->timer_3); // find neighbors time
	if(params->verbosity > 0)
		printf("Progress Indicator (no. of records): ");
	for(i=0; i < nrows; i++){
		// find matches for row i
		mmjFindMatches(i, params, cands, invIdx, lengths,
				docs, rsizes, rwgts, rfiwgts, cwgts, endptr, marker, marker2);

#ifdef EXTRATIMES
		gk_startwctimer(params->timer_4); // indexing time
		mmjIndexRow(i, params, docs, endptr, invIdx, rwgts, rfiwgts, cwgts);
		gk_stopwctimer(params->timer_4); // indexing time
#else
		mmjIndexRow(i, params, docs, endptr, invIdx, rwgts, rfiwgts, cwgts);
#endif

		if ( params->verbosity > 0 && i % progressInd == 0 ) {
			printf("%d%%..", (int) (100*i/(float)nrows));
			fflush(stdout);
		}
	}
	if(params->verbosity > 0)
		printf("100%%\n");
#ifdef EXTRACOUNTS
	for(j=0; j < ncols; j++)
		params->indexSize += invIdx->ends[j] - invIdx->starts[j];
#endif

	gk_stopwctimer(params->timer_3); // find neighbors time

	simSearchFinalize(params);
	gk_free((void**)&rsizes, &csizes, &rwgts, &rfiwgts, &cwgts, &cands, &endptr,
			&marker, &marker2, &lengths,
			&invIdx->ids, &invIdx->vals, &invIdx->ends, &invIdx->starts,
			&invIdx->accum, &invIdx->lens, &invIdx, LTERM);

}

/**
 * Index part of the vector after its search is complete
 */
void mmjIndexRow(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdxJ_t *invIdx,
		val_t *rwgts, val_t *rfiwgts, val_t *cwgts){
	idx_t j, k;
	val_t myval, maxrval;
	accum_t b1, b2, l;
	accum_t simT;

	if(params->sim == DA_SIM_TAN){
	    simT = 2*params->simT/(1 + params->simT);
	} else {
	    simT = params->simT;
	}

	ptr_t *rowptr = docs->rowptr; // index in rowind/rowval where each row starts
	idx_t *rowind = docs->rowind; // colid for each nnz
	val_t *rowval = docs->rowval; // val for each nnz

	for ( b1=0.0, j=rowptr[rid]; j < rowptr[rid+1]; j++ )
		b1 += rowval[j] * gk_min(rwgts[rid], cwgts[rowind[j]]);

	// index terms in row i
	for ( b2=1.0, l=0.5, j=rowptr[rid], endptr[rid]=rowptr[rid]; j < rowptr[rid+1]; j++ ) {
		myval = rowval[j];
		k = rowind[j];
		l -= 0.5*myval*myval;
		invIdx->ids[ invIdx->ends[k] ] = rid;
		invIdx->vals[ invIdx->ends[k] ] = rowval[j];
		invIdx->lens[ invIdx->ends[k] ] = l;
		invIdx->ends[k]++;
		endptr[rid]++;
		b1 -= myval * gk_min(rwgts[rid], cwgts[k]);
		b2 -= 0.5*myval*myval;
		if(gk_min(b1,b2) < simT )
			break;
	}

	// update max row val for inverse index
	for (maxrval = 0.0, j++; j < rowptr[rid+1]; j++ )
		if(rowval[j] > maxrval)
			maxrval = rowval[j];
	rfiwgts[rid] = maxrval;
}


/**
 * Process identified candidates
 */
idx_t mmjProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, val_t *lengths, ptr_t *endptr,
		idx_t ncands, idx_t *cands, accum_t *accum, val_t *marker, val_t *marker2, val_t *rwgts, val_t *rfiwgts){

	idx_t i, j, cid, dps = 0, um=0;
	accum_t simT, cval;
	ptr_t *rowptr;
	idx_t *rowind;
	val_t *rowval;
	da_csr_t *neighbors;

    if(params->sim == DA_SIM_TAN){
        simT = 2*((double)params->simT)/(1.0 + ((double)params->simT)) - TANEPS;
    } else {
        simT = params->simT;
    }

	rowptr = docs->rowptr; // index in rowind/rowval where each row starts
	rowind = docs->rowind; // colid for each nnz
	rowval = docs->rowval; // val for each nnz

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

	for ( i=0; i < ncands; i++ )
	{
		cid = cands[i];

		if ( accum[cid] + gk_min(rowptr[cid+1]-endptr[cid],
				rowptr[rid+1]-rowptr[rid]) *
				rfiwgts[cid] * rwgts[rid] >= simT )
		{
			um = 0;
			cval = accum[cid];
			for(j=endptr[cid]; j < rowptr[cid+1]; j++){
				if(marker[rowind[j]] > 0.0){
					cval += marker[rowind[j]] * rowval[j];
					if(um > 1){
						if(cval + marker2[rowind[j]] + lengths[j] < simT){
							cval = 0.0;
							break;
						}
						um = 0;
					}
				} else
					um++;
			}
			if(j == rowptr[cid+1])
				dps++;

			if(cval >= simT){
				if(params->sim == DA_SIM_TAN){
				    //convert to tanimoto from cosine and check initial threshod
				    cval /= ((docs->rnorms[rid] + docs->rnorms[cid])/
				            (sqrt(docs->rnorms[rid])*sqrt(docs->rnorms[cid])) - cval);
				    if(cval >= params->simT - TANEPS){
	                    //add to neighbors matrix
	                    params->nSimPairs++;
	                    da_storeSim(params, rid, cid, cval);
				    }
				} else {
                    //add to neighbors matrix
				    params->nSimPairs++;
                    da_storeSim(params, rid, cid, cval);
				}
			}
		}
		accum[cid] = -1;
	}

	return dps;
}




/**
 * Find matches for row id rid
 */
void mmjFindMatches(idx_t rid, params_t *params, idx_t *cands,
		da_invIdxJ_t *idx, val_t *lengths, da_csr_t *docs, idx_t *rsizes, val_t *rwgts, val_t *rfiwgts, val_t *cwgts,
		ptr_t *endptr, val_t *marker, val_t *marker2){

#ifdef EXTRATIMES
	gk_startwctimer(params->timer_1); // find candidates time
#endif
	ssize_t i, j;
	idx_t k, ncands, cid;
	ptr_t *rowptr, *starts, *ends;
	idx_t *rowind, *idxids;
	accum_t simT, v, myval, cval, add, upperbound, minsize;
	val_t *rowval, *idxvals, *idxlens;
	accum_t *accum, l, b1, b2;

	rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
	rowind  = docs->rowind; // col id for each nnz
	rowval  = docs->rowval; // val for each nnz
	idxids  = idx->ids;
	idxvals = idx->vals;
	idxlens = idx->lens;
	starts  = idx->starts;  // starting points for lists in the inv index
	ends    = idx->ends;    // ending points for lists in the inv index
	accum   = idx->accum;   // accumulator array
	ncands  = 0;            // number of candidates for this doc
	simT    = params->simT;

	// compute remscore and minsize
	for(b1=0.0, i=rowptr[rid]; i<rowptr[rid+1]; i++)
		b1 += rowval[i] * cwgts[ rowind[i] ]; // val * max col weight
	minsize = simT/rwgts[rid];

	for (b2=1, l=0.5, i=rowptr[rid]; i < rowptr[rid+1]; i++ )
	{
		k = rowind[i];
		myval = rowval[i];

		// remove those vectors from the inverted index that do
		// not have sufficient size.
		for (j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];
			if ( rsizes[cid] * rwgts[cid] < minsize ){
				starts[k]++;
			} else
				break;
		}

		l -= 0.5*myval*myval;

		for ( j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];  // potential candidate from the inv index

			if ( gk_min(b1,b2) >= simT || accum[cid] > 0.0 ){
				if ( accum[cid] == -1) {
					accum[cid] = 0.0;
					cands[ncands++] = cid;
				}
				accum[cid] +=  myval * idxvals[j];
				if(accum[cid] + l + idxlens[j] < simT)
					accum[cid] = 0.0;
			}
		}

		b1 -= myval * cwgts[k];
		b2 = l + 0.5;
		lengths[i] = l;

		marker[k] = myval;
		marker2[k] = l;
	}
#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_1); // find candidates time

	gk_startwctimer(params->timer_2); // process candidates time
#endif
	params->nCandidates  += ncands;
	params->nDotProducts += mmjProcessCandidates(params, rid, docs, lengths, endptr, ncands,
				cands, accum, marker, marker2, rwgts, rfiwgts);

#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_2); // process candidates time
	gk_startwctimer(params->timer_1); // find candidates time
#endif

	// reset marker
	for ( i=rowptr[rid]; i < rowptr[rid+1]; i++ ){
		k = rowind[i];
		marker[k] = marker2[k] = 0.0;
	}

#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_1); // find candidates time
#endif
}





/**
 * Reorder the document matrix in increasing number of column nnz
 * and decreasing order of row max val. Also compute and store
 * necessary stats about the data
 * Store initial order in rperm/cperm.
 */
void mmjReorderDocs(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
		val_t *cwgts, idx_t *rperm, idx_t *cperm){
	ssize_t i, j, k, nnz;
	idx_t nrows, ncols;
	da_iikv_t *orderColNnz = NULL;
	da_ivkv_t *orderMaxVal = NULL;
	char pv = 0;
	ptr_t *rowptr, *nrowptr;
	idx_t *rowind, *nrowind;
	val_t v, *rowval, *nrowval, *rnorms, *nrnorms=NULL;

	nrows  = (*docs)->nrows;
	ncols  = (*docs)->ncols;
	rowptr = (*docs)->rowptr;
	rowind = (*docs)->rowind;
	rowval = (*docs)->rowval;
	nnz    = rowptr[nrows];

	// for tan similarity, first compute and store norms, than normalize
	if(params->sim == DA_SIM_TAN){
	    da_csr_ComputeSquaredNorms(*docs, DA_ROW);
	    rnorms = (*docs)->rnorms;
	    // now normalize the rows
	    for(i=0; i < nrows; ++i){
	        v = 1.0/sqrt(rnorms[i]);
	        for(j=rowptr[i]; j < rowptr[i+1]; ++j){
	            rowval[j] *= v;
	        }
	    }
	    nrnorms = da_vmalloc(nrows, "mmjReorderDocs: nrnorms");
	}
    rnorms = (*docs)->rnorms;

	//record col nnzs
	da_iset(ncols, 0, csizes);
	for (i=0; i<nnz; i++)
		csizes[rowind[i]]++;

	//assign memory for ordering
	orderMaxVal = da_ivkvmalloc(nrows, "mmjReorderDocs: orderMaxVal");
	nrowptr = da_pmalloc(nrows+1, "mmjReorderDocs: nrowptr");
	nrowind = da_imalloc(nnz, "mmjReorderDocs: nrowind");
	nrowval = da_vmalloc(nnz, "mmjReorderDocs: nrowval");
	orderColNnz = da_iikvmalloc(ncols, "mmjReorderDocs: orderColNnz");

	//get new column order
	for(j=0; j < ncols; j++){
		orderColNnz[j].key = j;
		orderColNnz[j].val = csizes[j];
	}
	da_iikvsorti(ncols, orderColNnz); // sort columns
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
		orderMaxVal[i].key = i;
		orderMaxVal[i].val = rwgts[i];
	}
	da_ivkvsortd(nrows, orderMaxVal); // sort rows
	for(i=0; i < nrows; i++){
		rperm[i] = orderMaxVal[i].key; // store new order
		rwgts[i] = orderMaxVal[i].val; // new row weight after permuting
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
	if(params->sim == DA_SIM_TAN){
	    for(i=0; i < nrows; ++i){
	        nrnorms[i] = rnorms[rperm[i]];
	    }
	    gk_free((void**)&(*docs)->rnorms, LTERM);
	    (*docs)->rnorms = nrnorms;
	}
	(*docs)->rowptr = nrowptr;
	(*docs)->rowind = nrowind;
	(*docs)->rowval = nrowval;

	da_csr_SortIndices((*docs), DA_ROW);

	gk_free((void**)&orderColNnz, &orderMaxVal, &cperm, LTERM);

	params->rperm = rperm;

}
