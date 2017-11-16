/*!
 \file  ap.c
 \brief This file contains allPairs related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
void apFindMatches(idx_t rid, params_t *params, idx_t *candidates,
	da_invIdx_t *idx, da_csr_t *docs, idx_t *rsizes, val_t *rwgts, val_t *cwgts,
	ptr_t *endptr, val_t *marker);
void apIndexRow(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdx_t *invIdx,
		val_t *rwgts, val_t *cwgts);
idx_t apProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, ptr_t *endptr,
		idx_t cands, idx_t *candidates, accum_t *accum, val_t *marker, val_t *rfiwgts);
void apReorderDocs(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
	val_t *cwgts, idx_t *rperm, idx_t *cperm);
void ap2FindMatches(idx_t rid, params_t *params, idx_t *candidates,
		da_invIdx_t *idx, da_csr_t *docs, idx_t *rsizes, val_t *rwgts, val_t *cwgts,
		ptr_t *endptr, val_t *marker);
idx_t ap2ProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, ptr_t *endptr,
		idx_t cands, idx_t *candidates, accum_t *accum, val_t *marker, val_t *rfiwgts);



/**
 * Main entry point to AllPairs.
 */
void apFindNeighbors(params_t *params){

	ssize_t i, j, k, nnz;
	idx_t nrows, ncols, mrlen, mclen, progressInd;
	ptr_t *rowptr, *endptr;
	idx_t *candidates, *rsizes, *csizes, *rperm, *cperm;
	val_t b, maxrval, myval;
	val_t *rwgts = NULL, *cwgts = NULL, *marker;
	da_csr_t *docs = params->docs;
	da_invIdx_t *invIdx = NULL; //alternative inverted index for value-driven similarities//

	assert(docs->rowval);

	// allocate memory
	simSearchSetup(params);
	nrows  = docs->nrows;   // num rows
	ncols  = docs->ncols;   // num cols
	progressInd = ceil(nrows/(float)20);

	rsizes = da_imalloc(nrows, "apFindNeighbors: rsizes");
	csizes = da_ismalloc(ncols, 0.0, "apFindNeighbors: csizes");
	rperm  = da_imalloc(nrows, "apFindNeighbors: rperm");
	cperm  = da_imalloc(ncols, "apFindNeighbors: cperm");
	rwgts  = da_vsmalloc(nrows, 0.0, "apFindNeighbors: rwgts");
	cwgts  = da_vsmalloc(ncols, 0.0, "apFindNeighbors: cwgts");
	marker = da_vsmalloc(ncols, 0, "apFindNeighbors: marker");

	// reorder matrix
	apReorderDocs(params, &docs, rsizes, csizes, rwgts, cwgts, rperm, cperm);

	nnz    = docs->rowptr[nrows]; // nnz in docs
	rowptr = docs->rowptr; // index in rowind/rowval where each row starts

#ifdef EXTRACOUNTS
	// find the average length of rows & columns
	mrlen = ceil(da_isum(nrows, rsizes, 1) / (val_t) nrows);
	mclen = ceil(da_isum(ncols, csizes, 1) / (val_t) ncols);
	if(params->verbosity > 0)
		printf("mrlen: " PRNT_IDXTYPE ", mclen: " PRNT_IDXTYPE ".\n", mrlen, mclen);
#endif

	// allocate memory for candidates and the inverted indexes
	candidates = da_imalloc(nrows, "apFindNeighbors: candidates");
	endptr     = da_pmalloc(nrows, "apFindNeighbors: endptr");
	if(params->mode == MODE_AP){
		for(i=0; i < nrows; i++)
			endptr[i] = rowptr[i+1];
	}
	invIdx = gk_malloc(sizeof(da_invIdx_t), "apFindNeighbors: invIdx");
	invIdx->ids    = da_ismalloc(nnz, 0, "apFindNeighbors: invIdx->ids");
	invIdx->vals   = da_vsmalloc(nnz, 0.0, "apFindNeighbors: invIdx->ids");
	invIdx->starts = da_pmalloc(ncols, "apFindNeighbors: invIdx->starts");
	invIdx->ends   = da_pmalloc(ncols, "apFindNeighbors: invIdx->ends");
	for(j=0, k=0; j < ncols; j++){
		invIdx->starts[j] = invIdx->ends[j] = k;
		k += csizes[j];
	}
	invIdx->accum = gk_malloc(nrows*sizeof(accum_t), "apFindNeighbors: invIdx->accum");
	memset(invIdx->accum, 0, nrows*sizeof(accum_t));

	// for each x \in V do
	gk_startwctimer(params->timer_3); // find neighbors time
	if(params->verbosity > 0)
		printf("Progress Indicator (no. of records): ");
	if(params->mode == MODE_AP2)
		for(i=0; i < nrows; i++){
			// find matches for row i
			ap2FindMatches(i, params, candidates, invIdx,
					docs, rsizes, rwgts, cwgts, endptr, marker);

	#ifdef EXTRATIMES
			gk_startwctimer(params->timer_4); // indexing time
			apIndexRow(i, params, docs, endptr, invIdx, rwgts, cwgts);
			gk_stopwctimer(params->timer_4); // indexing time
	#else
			apIndexRow(i, params, docs, endptr, invIdx, rwgts, cwgts);
	#endif

			if ( params->verbosity > 0 && i % progressInd == 0 ) {
				printf("%d%%..", (int) (100*i/(float)nrows));
				fflush(stdout);
			}
		}
	else
		for(i=0; i < nrows; i++){
			// find matches for row i
			apFindMatches(i, params, candidates, invIdx,
					docs, rsizes, rwgts, cwgts, endptr, marker);

	#ifdef EXTRATIMES
			gk_startwctimer(params->timer_4); // indexing time
			apIndexRow(i, params, docs, endptr, invIdx, rwgts, cwgts);
			gk_stopwctimer(params->timer_4); // indexing time
	#else
			apIndexRow(i, params, docs, endptr, invIdx, rwgts, cwgts);
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
	params->nPruneDotP = params->nCandidates - params->nDotProducts;
#endif

	simSearchFinalize(params);
	gk_stopwctimer(params->timer_3); // find neighbors time

	gk_free((void**)&rsizes, &csizes, &rwgts, &cwgts, &candidates, &endptr, &marker,
			&invIdx->ids, &invIdx->vals, &invIdx->ends, &invIdx->starts,
			&invIdx->accum, &invIdx, LTERM);

}

/**
 * Index part of the vector after its search is complete
 */
void apIndexRow(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdx_t *invIdx,
		val_t *rwgts, val_t *cwgts){
	idx_t j;
	val_t myval, maxrval;
	double b;
	ptr_t *rowptr = docs->rowptr; // index in rowind/rowval where each row starts
	idx_t *rowind = docs->rowind; // colid for each nnz
	val_t *rowval = docs->rowval; // val for each nnz

	if(params->mode == MODE_AP2){

		for ( b=0.0, j=rowptr[rid]; j < rowptr[rid+1]; j++ )
			b += rowval[j] * gk_min(rwgts[rid], cwgts[rowind[j]]);

		// index terms in row i
		for ( j=rowptr[rid], endptr[rid]=rowptr[rid]; j < rowptr[rid+1]; j++ ) {
			myval = rowval[j];
			b -= (double) myval * gk_min(rwgts[rid], cwgts[rowind[j]]);
			invIdx->ids[ invIdx->ends[rowind[j]] ] = rid;
			invIdx->vals[ invIdx->ends[rowind[j]] ] = rowval[j];
			invIdx->ends[ rowind[j] ]++;

			if(b < params->simT )
				break;
		}
		endptr[rid] = ++j;

		// update max row val for forward index
		for (maxrval = 0.0; j < rowptr[rid+1]; j++ )
			if(rowval[j] > maxrval)
				maxrval = rowval[j];
		rwgts[rid] = maxrval;

		return;
	}

	// index terms in row i
	maxrval = b = 0.0;
	for ( j = rowptr[rid]; j < endptr[rid]; j++ ) {
		myval = rowval[j];
		b += myval * gk_min(rwgts[rid], cwgts[rowind[j]]);
		if(b >= params->simT )
			break;
		if( myval > maxrval )
			maxrval = myval;
	}
	// truncate the rest of the vector, since we're going to
	// index it.
	endptr[rid] = j;
	// update max row val for row rid
	rwgts[rid] = maxrval;
	// index the remaining part of the vector
	for ( ; j < rowptr[rid+1]; j++ ){
		invIdx->ids[ invIdx->ends[rowind[j]] ] = rid;
		invIdx->vals[ invIdx->ends[rowind[j]] ] = rowval[j];
		invIdx->ends[rowind[j]]++;
	}
}


/**
 * Process identified candidates
 */
idx_t apProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, ptr_t *endptr,
		idx_t ncands, idx_t *candidates, accum_t *accum, val_t *marker, val_t *rwgts){

	idx_t i, j, cid, dps = 0;
	val_t upperbound;
	accum_t cval;
	ptr_t *rowptr;
	idx_t *rowind;
	val_t *rowval;
	da_csr_t *neighbors;

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
		cid = candidates[i];

		upperbound = gk_min(endptr[cid]-rowptr[cid],
				endptr[rid]-rowptr[rid]) *
				rwgts[cid] * rwgts[rid] + accum[cid];
		if ( upperbound < params->simT )
		{
			accum[cid] = 0;
			continue;
		}

		// partial dot product between rid and cid
		cval = accum[cid];
		for(j=rowptr[cid]; j < endptr[cid]; j++)
				cval += marker[rowind[j]] * rowval[j];
		dps++;

		if ( cval >= params->simT ){
			params->nSimPairs++;
			//add to neighbors matrix
			da_storeSim(params, rid, cid, cval);
		}
		accum[cid] = 0;
	}

	return dps;
}


/**
 * Process identified candidates - alternate order of cols
 */
idx_t ap2ProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, ptr_t *endptr,
		idx_t ncands, idx_t *candidates, accum_t *accum, val_t *marker, val_t *rwgts){

	idx_t i, j, cid, dps = 0;
	val_t upperbound;
	accum_t cval;
	ptr_t *rowptr;
	idx_t *rowind;
	val_t *rowval;
	da_csr_t *neighbors;

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
		cid = candidates[i];

		upperbound = gk_min(rowptr[cid+1]-endptr[cid],
				rowptr[rid+1]-rowptr[rid]) *
				rwgts[cid] * rwgts[rid] + accum[cid];
		if ( upperbound < params->simT )
		{
			accum[cid] = 0;
			continue;
		}

		// partial dot product between rid and cid
		cval = accum[cid];
		for(j=endptr[cid]; j < rowptr[cid+1]; j++){
			cval += (double) marker[rowind[j]] * rowval[j];
		}
		dps++;

		if ( cval >= params->simT ){
			params->nSimPairs++;
			//add to neighbors matrix
			da_storeSim(params, rid, cid, cval);
		}
		accum[cid] = 0;
	}

	return dps;
}


/**
 * Find matches for row id rid
 */
void apFindMatches(idx_t rid, params_t *params, idx_t *candidates,
		da_invIdx_t *idx, da_csr_t *docs, idx_t *rsizes, val_t *rwgts, val_t *cwgts,
		ptr_t *endptr, val_t *marker){

#ifdef EXTRATIMES
	gk_startwctimer(params->timer_1); // find candidates time
#endif
	ssize_t i, j, b;
	idx_t k, minsize, cands, ncols, cid;
	ptr_t *rowptr, *starts, *ends;
	idx_t *rowind, *idxids;
	val_t simT, v, remscore, myval, cval, add, upperbound;
	val_t *rowval, *idxvals;
	accum_t *accum;

	ncols   = docs->ncols;  // num cols
	rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
	rowind  = docs->rowind; // col id for each nnz
	rowval  = docs->rowval; // val for each nnz
	idxids  = idx->ids;
	idxvals = idx->vals;
	starts  = idx->starts;  // starting points for lists in the inv index
	ends    = idx->ends;    // ending points for lists in the inv index
	accum   = idx->accum;   // accumulator array
	cands   = 0;            // number of candidates for this doc
	simT    = params->simT;

	// compute remscore and minsize
	for(i=rowptr[rid], remscore = 0.0; i<rowptr[rid+1]; i++)
		remscore += rowval[i] * cwgts[ rowind[i] ]; // val * max col weight
	minsize = ceil(simT/rwgts[rid]);

	if(remscore < simT || minsize > ncols){
#ifdef EXTRATIMES
		gk_stopwctimer(params->timer_1); // find candidates time
#endif
		return;
	}


	for ( i=rowptr[rid+1]-1; i >= rowptr[rid]; i-- )
	{
		k = rowind[i];
		myval = rowval[i];
		if(i < endptr[rid])
			marker[k] = myval;

		// remove those vectors from the inverted index that do
		// not have sufficient size.
		for (b = starts[k], j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];
			if ( rsizes[cid] >= minsize )
				break;
		}
		starts[k] += j-b;
#ifdef EXTRACOUNTS
		params->nPruneMinsize += j-b;
#endif

		for ( j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];  // potential candidate from the inv index

			if ( remscore >= simT || accum[cid] > 0 ){
				add = myval * idxvals[j];
				if ( accum[cid] <= 0 && add > 0) {
					candidates[cands++] = cid;
				}
				accum[cid] += add;
			}
		}

		remscore -= myval * cwgts[k];
	}
#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_1); // find candidates time

	gk_startwctimer(params->timer_2); // process candidates time
#endif
	params->nCandidates  += cands;
	params->nDotProducts += apProcessCandidates(params, rid, docs, endptr, cands,
				candidates, accum, marker, rwgts);
#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_2); // process candidates time
#endif

#ifdef EXTRATIMES
	gk_startwctimer(params->timer_1); // find candidates time
#endif
	// reset marker
	for ( i=rowptr[rid]; i < endptr[rid]; i++ )
		marker[rowind[i]] = 0;

#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_1); // find candidates time
#endif
}



/**
 * Find matches for row id rid - alternate column order
 */
void ap2FindMatches(idx_t rid, params_t *params, idx_t *candidates,
		da_invIdx_t *idx, da_csr_t *docs, idx_t *rsizes, val_t *rwgts, val_t *cwgts,
		ptr_t *endptr, val_t *marker){

#ifdef EXTRATIMES
	gk_startwctimer(params->timer_1); // find candidates time
#endif
	ssize_t i, j, b;
	idx_t k, minsize, cands, ncols, cid;
	ptr_t *rowptr, *starts, *ends;
	idx_t *rowind, *idxids;
	val_t simT, v, remscore, myval, cval, add, upperbound;
	val_t *rowval, *idxvals;
	accum_t *accum;

	ncols   = docs->ncols;  // num cols
	rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
	rowind  = docs->rowind; // col id for each nnz
	rowval  = docs->rowval; // val for each nnz
	idxids  = idx->ids;
	idxvals = idx->vals;
	starts  = idx->starts;  // starting points for lists in the inv index
	ends    = idx->ends;    // ending points for lists in the inv index
	accum   = idx->accum;   // accumulator array
	cands   = 0;            // number of candidates for this doc
	simT    = params->simT;

	// compute remscore and minsize
	for(i=rowptr[rid], remscore = 0.0; i<rowptr[rid+1]; i++)
		remscore += rowval[i] * cwgts[ rowind[i] ]; // val * max col weight
	minsize = ceil(simT/rwgts[rid]);

	if(remscore < simT || minsize > ncols){
#ifdef EXTRATIMES
		gk_stopwctimer(params->timer_1); // find candidates time
#endif
		return;
	}


	for ( i=rowptr[rid]; i < rowptr[rid+1]; i++ )
	{
		k = rowind[i];
		myval = rowval[i];
		marker[k] = myval;

		// remove those vectors from the inverted index that do
		// not have sufficient size.
		for (b = starts[k], j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];
			if ( rsizes[cid] >= minsize )
				break;
		}
		starts[k] += j-b;
#ifdef EXTRACOUNTS
		params->nPruneMinsize += j-b;
#endif

		for ( j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];  // potential candidate from the inv index

			if ( remscore >= simT || accum[cid] > 0 ){
				add = myval * idxvals[j];
				if ( accum[cid] <= 0 && add > 0) {
					candidates[cands++] = cid;
				}
				accum[cid] += add;
			}
		}

		remscore -= myval * cwgts[k];
	}
#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_1); // find candidates time

	gk_startwctimer(params->timer_2); // process candidates time
#endif
	params->nCandidates  += cands;
	params->nDotProducts += ap2ProcessCandidates(params, rid, docs, endptr, cands,
				candidates, accum, marker, rwgts);
#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_2); // process candidates time

	gk_startwctimer(params->timer_1); // find candidates time
#endif
	// reset marker
	for ( i=rowptr[rid]; i < rowptr[rid+1]; i++ )
		marker[rowind[i]] = 0;

#ifdef EXTRATIMES
	gk_stopwctimer(params->timer_1); // find candidates time
#endif
}


/**
 * Reorder the document matrix in decreasing number of column nnz
 * and decreasing(ap)/increasing(ap2) order of row max val. Also compute and store
 * necessary stats about the data
 * Store initial order in rperm/cperm.
 */
void apReorderDocs(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
		val_t *cwgts, idx_t *rperm, idx_t *cperm){
	ssize_t i, j, k, nnz;
	idx_t nrows, ncols;
	da_iikv_t *orderColNnz = NULL;
	da_ivkv_t *orderMaxVal = NULL;
	char pv = 0;
	ptr_t *rowptr, *nrowptr;
	idx_t *rowind, *nrowind, *tperm;
	val_t v, *rowval, *nrowval;

	nrows  = (*docs)->nrows;
	ncols  = (*docs)->ncols;
	rowptr = (*docs)->rowptr;
	rowind = (*docs)->rowind;
	rowval = (*docs)->rowval;
	nnz    = rowptr[nrows];

	//record col nnzs
	da_iset(ncols, 0, csizes);
	for (i=0; i<nnz; i++)
		csizes[rowind[i]]++;

	//assign memory for ordering
	orderMaxVal = da_ivkvmalloc(nrows, "apReorderDocs: orderMaxVal");
	nrowptr = da_pmalloc(nrows+1, "apReorderDocs: nrowptr");
	nrowind = da_imalloc(nnz, "apReorderDocs: nrowind");
	nrowval = da_vmalloc(nnz, "apReorderDocs: nrowval");

	orderColNnz = da_iikvmalloc(ncols, "apReorderDocs: orderColNnz");
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
	(*docs)->rowptr = nrowptr;
	(*docs)->rowind = nrowind;
	(*docs)->rowval = nrowval;

	da_csr_SortIndices((*docs), DA_ROW);

	//inverse the row permutation
	tperm = rperm;
	rperm = da_imalloc(nrows, "l2apReorderDocs: rperm");
	for ( i=0; i<nrows; i++ )
		rperm[i] = tperm[i];

	gk_free((void**)&orderColNnz, &orderMaxVal, &cperm, &tperm, LTERM);

	params->rperm = rperm;

}
