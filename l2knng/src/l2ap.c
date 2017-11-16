/*!
 \file  l2ap.c
 \brief This file contains l2ap related functions

 \author David C. Anastasiu
 */

#include "includes.h"

// forward declarations
idx_t l2apProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, val_t *lengths, val_t *sums,
		ptr_t *endptr, idx_t ncands, idx_t *cands, accum_t *accum,
		val_t *hashval, val_t *hashwgt, val_t *hashlen, val_t *hashsum, idx_t *hashsz, val_t *rwgts, val_t *rfiwgts,
		val_t *ps);
void l2apIndexRow(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdxJ_t *invIdx,
		val_t *rwgts, val_t *rfiwgts, val_t *cwgts, val_t *criwgts, val_t *hashwgt, val_t *ps);
void l2apFindMatches(idx_t rid, params_t *params, idx_t *cands,
		da_invIdxJ_t *idx, val_t *lengths, val_t *sums, da_csr_t *docs, idx_t *rsizes,
		val_t *rwgts, val_t *rfiwgts, val_t *cwgts, val_t *criwgts, val_t *ps, ptr_t *endptr,
		val_t *hashval, val_t *hashwgt, val_t *hashlen, val_t *hashsum, idx_t *hashsz);

/**
 * Main entry point to L2AP.
 */
void l2apFindNeighbors(params_t *params){

	ssize_t i, j, k, nnz;
	idx_t nrows, ncols, mrlen, mclen, pct, progressInd;
	ptr_t *rowptr, *endptr;
	idx_t *cands, *rsizes, *csizes, *rperm, *cperm, *hashsz = NULL;
	val_t maxrval, myval;
	val_t *rwgts = NULL, *rfiwgts = NULL, *hashwgt = NULL, *cwgts = NULL,
			*criwgts = NULL, *hashval = NULL, *hashlen = NULL, *hashsum = NULL,
			*ps = NULL, *lengths = NULL, *sums = NULL;
	da_invIdxJ_t *invIdx = NULL; //alternative inverted index for value-driven similarities//
	da_csr_t *docs = params->docs;
	char *querydocs = params->querydocs;

	nrows  = docs->nrows;   // num rows
	ncols  = docs->ncols;   // num cols
	params->simT2 = params->simT * params->simT;
	progressInd = ceil(nrows/(float)20);

    rsizes  = params->rsizes;
    if(!rsizes)
        rsizes = params->rsizes = da_imalloc(nrows, "l2apFindNeighbors: rsizes");
	csizes  = params->csizes;
    if(!csizes)
        csizes = params->csizes = da_ismalloc(ncols, 0, "l2apFindNeighbors: csizes");
	rperm   = params->rperm;
    if(!rperm)
        rperm = params->rperm = da_imalloc(nrows, "l2apFindNeighbors: rperm");
	cperm   = params->cperm;
    if(!cperm)
        cperm = params->cperm = da_imalloc(ncols, "l2apFindNeighbors: cperm");
	rwgts   = params->rwgts;
    if(!rwgts)
        rwgts = params->rwgts = da_vsmalloc(nrows, 0.0, "l2apFindNeighbors: rwgts");
	rfiwgts = da_vmalloc(nrows, "l2apFindNeighbors: rfiwgts");
	cwgts   = params->cwgts;
    if(!cwgts)
        cwgts = params->cwgts = da_vsmalloc(ncols, 0.0, "l2apFindNeighbors: cwgts");
#if defined(RS1) && defined(RS3)
	criwgts = da_vsmalloc(ncols, 0.0, "l2apFindNeighbors: criwgts");
#endif
	hashval = da_vsmalloc(ncols, 0, "l2apFindNeighbors: hashval");      // hash values of x
#if defined(DP7) || defined(DP8)
	hashsz  = da_imalloc(ncols, "l2apFindNeighbors: hashsz");
#endif
#if defined(L2CV) || defined(LENCV)
	hashlen = da_vsmalloc(ncols, 0, "l2apFindNeighbors: hashlen");      // hash suffix lengths of x
#endif
#if defined(DP3) || defined(DP4)
	hashsum = da_vsmalloc(ncols, 0, "l2apFindNeighbors: hashsum");      // hash suffix sums of x
#endif
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
	hashwgt = da_vmalloc(ncols, "l2apFindNeighbors: hashwgt");
#endif
#ifdef PSCV
	ps      = da_vmalloc(nrows, "l2ap2FindNeighbors: ps");        // score prior to indexing threshold for each row
#endif
	if(querydocs == NULL)
	    querydocs = params->querydocs = da_csmalloc(nrows, 1, "l2ap2FindNeighbors: querydocs");

	/** matrix already reordered for kl2ap job */
/**
 *  if(preprocess)
 *	    l2apReorderDocs(params, &docs, rsizes, csizes, rwgts, cwgts, rperm, cperm);
**/

	nnz    = docs->rowptr[nrows]; // nnz in docs
	rowptr = docs->rowptr; // index in rowind/rowval where each row starts

	// allocate memory for candidates and the inverted indexes -- needs nnz & col counts
	cands  = da_imalloc(nrows, "l2apFindNeighbors: candidates");
#if defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4)
	sums    = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: sums");      // keep track of suffix sums
#endif
	lengths = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: lengths");   // keep track of suffix lengths
	endptr = da_pmalloc(nrows, "l2apFindNeighbors: endptr");
	invIdx = da_malloc(sizeof(da_invIdxJ_t), "l2apFindNeighbors: invIdx");
	invIdx->ids    = da_ismalloc(nnz, 0, "l2apFindNeighbors: invIdx->ids");
	invIdx->vals   = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: invIdx->ids");
	invIdx->lens   = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: invIdx->ids");
	invIdx->starts = da_pmalloc(ncols, "l2apFindNeighbors: invIdx->starts");
	invIdx->ends   = da_pmalloc(ncols, "l2apFindNeighbors: invIdx->ends");
	for(j=0, k=0; j < ncols; j++){
		invIdx->starts[j] = invIdx->ends[j] = k;
		k += csizes[j];
	}
	invIdx->accum = da_malloc(nrows*sizeof(accum_t), "l2apFindNeighbors: invIdx->accum");
	for(i=0; i < nrows; i++)
		invIdx->accum[i] = -1;

	// for each x \in V do
	da_startwctimer(params->timer_3); // find neighbors time
	if(params->verbosity > 0)
		printf("Searching... ");
    da_progress_init(pct, progressInd, nrows);
	for(i=0; i < nrows; i++){
		endptr[i] = rowptr[i+1];
		// find matches for row i
		l2apFindMatches(i, params, cands, invIdx, lengths, sums,
				docs, rsizes, rwgts, rfiwgts, cwgts, criwgts, ps, endptr,
				hashval, hashwgt, hashlen, hashsum, hashsz);

#ifdef EXTRATIMES
		da_startwctimer(params->timer_4); // indexing time
		l2apIndexRow(i, params, docs, endptr, invIdx, rwgts, rfiwgts, cwgts, criwgts, hashwgt, ps);
		da_stopwctimer(params->timer_4); // indexing time
#else
		l2apIndexRow(i, params, docs, endptr, invIdx, rwgts, rfiwgts, cwgts, criwgts, hashwgt, ps);
#endif

		if ( params->verbosity > 0 && i % progressInd == 0 )
		    da_progress_advance(pct);

	}
    if (params->verbosity > 0)
        da_progress_finalize(pct);
#ifdef EXTRACOUNTS
	for(j=0; j < ncols; j++)
		params->indexSize += invIdx->ends[j] - invIdx->starts[j];
#endif
	da_stopwctimer(params->timer_3); // find neighbors time

/**	simSearchFinalize(params); **/

	da_free((void**)&rfiwgts, &criwgts,
			&cands, &endptr, &hashval, &hashwgt, &hashlen, &hashsum, &hashsz,
			&lengths, &sums, &ps,
			&invIdx->ids, &invIdx->vals, &invIdx->ends, &invIdx->starts,
			&invIdx->accum, &invIdx->lens, &invIdx, LTERM);

}

/**
 * Index part of the vector after its search is complete
 */
void l2apIndexRow(idx_t rid, params_t *params, da_csr_t *docs, ptr_t *endptr, da_invIdxJ_t *invIdx,
		val_t *rwgts, val_t *rfiwgts, val_t *cwgts, val_t *criwgts, val_t *hashwgt, val_t *ps){
	ssize_t j, k;
	val_t myval, maxrval;
	double b1, b2, b2l, sqb2, pscore, simT;
	ptr_t *rowptr = docs->rowptr; // index in rowind/rowval where each row starts
	idx_t *rowind = docs->rowind; // colid for each nnz
	val_t *rowval = docs->rowval; // val for each nnz
	simT = params->simT;
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
	for ( endptr[rid] = rowptr[rid+1], j = rowptr[rid]; j < endptr[rid]; j++ ) {
		myval = rowval[j];
		b1 += (double) myval * da_min(rwgts[rid], cwgts[rowind[j]]);
		// adjust b2 estimate of prefix similarity with all other vectors
#if defined(L2PS) || defined (IDXL2)
		b2 += (double) myval * myval;
		sqb2 = sqrt(b2);
#elif defined(LENPS) || defined (IDXLEN)
		b2 += (double) 0.5 * myval * myval;
#endif

		// check whether to start indexing
#ifdef L2PS
		if(da_min(b1,sqb2) >= simT )
			break;
#elif defined(LENPS)
		if(da_min(b1,b2) >= simT )
			break;
#else
		if(b1 >= simT )
			break;
#endif

#ifdef PSCV
		//update pscore
	#ifdef L2PS
		pscore = da_min(b1,sqb2);
	#elif defined(LENPS)
		pscore = da_min(b1,b2);
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

#if defined(RS1) && defined(RS3)
		// adjust inverted index max col weights
		if(myval > criwgts[rowind[j]])
			criwgts[rowind[j]] = myval;
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
#elif defined(LENPS)
	b2 -= 0.5;
#endif

	// index the remaining part of the vector
	for ( ; j < rowptr[rid+1]; j++ ){
		k = rowind[j];
#if defined (IDXL2)
		invIdx->lens[ invIdx->ends[k] ] = sqb2;
		b2 += (double) rowval[j] * rowval[j];
		sqb2 = sqrt(b2);
#elif defined(IDXLEN)
		invIdx->lens[ invIdx->ends[k] ] = b2;
		b2 += (double) 0.5 * rowval[j] * rowval[j];
#endif
		invIdx->ids[ invIdx->ends[k] ] = rid;
		invIdx->vals[ invIdx->ends[k] ] = rowval[j];
		invIdx->ends[k]++;
#if defined(RS1) && defined(RS3)
		if(rowval[j] > criwgts[k])
			criwgts[k] = rowval[j];
#endif
	}
}



idx_t l2apProcessCandidates(params_t *params, idx_t rid, da_csr_t *docs, val_t *lengths, val_t *sums,
		ptr_t *endptr, idx_t ncands, idx_t *cands, accum_t *accum,
		val_t *hashval, val_t *hashwgt, val_t *hashlen, val_t *hashsum, idx_t *hashsz, val_t *rwgts, val_t *rfiwgts,
		val_t *ps){

	ssize_t i, j, k, p, r, ls, sz, cid, nneighbs = -1, dps = 0;
	val_t rwgt, rsum;
	double cval, simT, simTe;
	ptr_t *rowptr;
	idx_t *rowind;
	val_t *rowval;

	rowptr = docs->rowptr; // index in rowind/rowval where each row starts
	rowind = docs->rowind; // colid for each nnz
	rowval = docs->rowval; // val for each nnz

	simT  = params->simT;
#if defined(L2CV)
	simTe = params->simT-LENBNDEPS;
#endif
#if defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7)
	rwgt  = rwgts[rid];
#endif
#if defined(DP1) || defined(DP2)
	rsum  = sums[rowptr[rid+1]-1];
#endif
#if defined(DP5) || defined(DP6)
	sz    = rowptr[rid+1]-rowptr[rid];
#endif

	for ( i=0; i < ncands; i++ )
	{
		cid = cands[i];

		if(accum[cid] <= 0.0){
			accum[cid] = -1;
			continue;
		}

		cval = accum[cid];
		accum[cid] = -1;

#ifdef PSCV
		if ( cval + ((double)ps[cid]) < simT ){
			#ifdef EXTRACOUNTS
				params->nPrunePscore++;
			#endif
			continue;
		}
#endif

#ifdef DP1
		if ( cval + ((double) da_min(rfiwgts[cid] * rsum,  rwgt * sums[endptr[cid]]) ) < simT ){
			#ifdef EXTRACOUNTS
				params->nPruneDotP++;
			#endif
			continue;
		}
#endif

#ifdef DP5
		if( cval + ((double) da_min(endptr[cid]-rowptr[cid], sz) * rfiwgts[cid] * rwgt) < simT){
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
		for(j=endptr[cid]-1, p=endptr[cid]-r; j >= r; j--, p--)
			if(hashval[rowind[j]] > 0.0)
				break;

		if(j < r)
			goto check;

	#ifdef DP2
		if ( cval + ((double) da_min(rfiwgts[cid] * rsum,  hashwgt[rowind[j]] * sums[j]) ) < simT ){
			#ifdef EXTRACOUNTS
				params->nPruneDotP2++;
			#endif
			continue;
		}
	#elif defined(DP3)
		if ( cval + ((double) da_min(rfiwgts[cid] * hashsum[rowind[j]],  rwgt * sums[j]) ) < simT ){
			#ifdef EXTRACOUNTS
				params->nPruneDotP2++;
			#endif
			continue;
		}
	#elif defined(DP4)
		if ( cval + ((double) da_min(rfiwgts[cid] * hashsum[rowind[j]],  hashwgt[rowind[j]] * sums[j]) ) < simT ){
			#ifdef EXTRACOUNTS
				params->nPruneDotP2++;
			#endif
			continue;
		}
	#elif defined(DP6)
		if ( cval + ((double) da_min(p, sz) * rfiwgts[cid] * hashwgt[rowind[j]]) < simT ){
			#ifdef EXTRACOUNTS
				params->nPruneDotP2++;
			#endif
			continue;
		}
	#elif defined(DP7)
		if ( cval + ((double) da_min(p, hashsz[rowind[j]]) * rfiwgts[cid] * rwgt) < simT ){
			#ifdef EXTRACOUNTS
				params->nPruneDotP2++;
			#endif
			continue;
		}
	#elif defined(DP8)
		if ( cval + ((double) da_min(p, hashsz[rowind[j]]) * rfiwgts[cid] * hashwgt[rowind[j]]) < simT ){
			#ifdef EXTRACOUNTS
				params->nPruneDotP2++;
			#endif
			continue;
		}
	#endif

#else
		j=endptr[cid]-1;
#endif

#ifdef L2CV
		for( ; j >= r; j--)
			if(hashval[rowind[j]] > 0.0){
				cval += (double) hashval[rowind[j]] * rowval[j];
				if(j > r && cval + ((double) hashlen[rowind[j]] * lengths[j]) < simTe){
					#ifdef EXTRACOUNTS
						params->nPruneLength2++;
					#endif
					goto nexcid;
				}
			}
#elif defined(LENCV)
		for( ; j >= r; j--)
			if(hashval[rowind[j]] > 0.0){
				cval += (double) hashval[rowind[j]] * rowval[j];
				if(j > r && cval + ((double) hashlen[rowind[j]]) + ((double) lengths[j]) < simT){
					#ifdef EXTRACOUNTS
						params->nPruneLength2++;
					#endif
					goto nexcid;
				}
			}
#else
		// partial dot product between rid and cid
		for( ; j >= r; j--)
			if(hashval[rowind[j]] > 0.0)
				cval += hashval[rowind[j]] * rowval[j];
#endif
		dps++;

		check:  // check similarity value

		if(cval >= simT){
			params->nSimPairs++;
			//add to KNNG
			da_storeSimKnn(params, rid, cid, cval);
		}

		nexcid:
			continue;

	}

	return dps;
}



/**
 * Find matches for row id rid
 */
void l2apFindMatches(idx_t rid, params_t *params, idx_t *cands,
		da_invIdxJ_t *idx, val_t *lengths, val_t *sums, da_csr_t *docs, idx_t *rsizes,
		val_t *rwgts, val_t *rfiwgts, val_t *cwgts, val_t *criwgts, val_t *ps, ptr_t *endptr,
		val_t *hashval, val_t *hashwgt, val_t *hashlen, val_t *hashsum, idx_t *hashsz ){

#ifdef EXTRATIMES
	da_startwctimer(params->timer_1); // find candidates time
#endif
	ssize_t i, j, s, ls;
	idx_t k, ncands, cid;
	ptr_t *rowptr, *starts, *ends;
	idx_t *rowind, *idxids;
	val_t v, minsize, myval;
	val_t *rowval, *idxvals, *idxlens;
	accum_t *accum;
	double simT, simTe, cval, sum, bnd, rs1, rs2, bsq, acc;
	char *qdocs;
	char qr;

	rowptr  = docs->rowptr; // index in rowind/rowval where each row starts
	rowind  = docs->rowind; // col id for each nnz
	rowval  = docs->rowval; // val for each nnz
	idxids  = idx->ids;
	idxvals = idx->vals;
	starts  = idx->starts;  // starting points for lists in the inv index
	ends    = idx->ends;    // ending points for lists in the inv index
	accum   = idx->accum;   // accumulator array
	qdocs   = params->querydocs; // docs we're interested in
	qr      = qdocs[rid];
	ncands  = 0;            // number of candidates for this doc
	simT    = params->simT;
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
	simTe   = params->simT-LENBNDEPS;
#endif
#ifdef RS1
	rs1 = 0.0;
#endif

	// define vars for hashing DP data
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
	v   = 0.0;
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
	for(i=rowptr[rid]; i<rowptr[rid+1]; i++){
		#ifdef RS1
			#ifdef RS3
			rs1 += rowval[i] * criwgts[ rowind[i] ]; // val * max col weight
			#else
			rs1 += rowval[i] * cwgts[ rowind[i] ]; // val * max col weight
			#endif
		#endif
		#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
			if(rowval[i] > v)
				v = rowval[i];
			hashwgt[rowind[i]] = v;  // prefix scan of max row values, i.e. remaining max weight in row
		#endif
		#if defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4)
			sum += rowval[i];
			sums[i] = sum;  // prefix sum
		#endif
		#if defined(DP3) || defined(DP4)
			hashsum[rowind[i]] = sums[i];  // hash sum
		#endif
		#if defined(DP7) || defined(DP8)
			hashsz[rowind[i]]  = k++;
		#endif
	}
	#if defined(RS1) && !defined(RS3)
	if(rs1 < simT){
		#ifdef EXTRATIMES
		da_stopwctimer(params->timer_1); // find candidates time
		#endif
		return;
	}
	#endif
#endif

	// find minimum size
#ifdef SZ3
	minsize = simT/rwgts[rid];
#endif

	for ( i=rowptr[rid+1]-1; i >= rowptr[rid]; i-- )
	{
		k = rowind[i];
		myval = rowval[i];
		hashval[k] = myval;

#ifdef SZ3
		// remove those vectors from the inverted index that do
		// not have sufficient size.
		for (s = starts[k], j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];
			if ( rsizes[cid]
					#ifndef SZ1
			            * rwgts[cid]
					#endif
			                    >= minsize )
				break;
		}
		starts[k] += j-s;
		#ifdef EXTRACOUNTS
		params->nPruneMinsize += j-s;
		#endif
#endif

#ifdef RS1
	#if defined(RS4)
		bnd   = da_min(rs1, bsq);
	#elif defined(RS2)
		bnd   = da_min(rs1, rs2);
	#else
		bnd = rs1;
	#endif
#else
	#if defined(RS4)
		bnd   = bsq;
	#elif defined(RS2)
		bnd   = rs2;
	#else
		bnd = 0;
		da_errexit("Invalid pruning selected at compile time. One of RS1-RS4 must be selected.");
	#endif
#endif

	#if defined(IDXL2)
		rs2  -= (double) myval * myval;
		bsq   = sqrt(rs2);
	#elif defined(IDXLEN)
		rs2  -= (double) 0.5 * myval * myval;
	#endif

	#if defined(L2CV)
		hashlen[k] = lengths[i] = bsq;
	#elif defined(LENCV)
		hashlen[k] = lengths[i] = rs2;
	#endif

		for ( j = starts[k]; j < ends[k]; j++ ){
			cid = idxids[j];  // potential candidate from the inv index
			acc = accum[cid];

			if ( acc > 0.0 ){
				// accumulate
				acc = accum[cid] += myval * idxvals[j];
				// check l2-norm bound
				#if defined(L2CG)
				if(acc + bsq * idxlens[j] < simTe){
					accum[cid] = -2;
					#ifdef EXTRACOUNTS
					params->nPruneLength++;
					#endif
				}
				#elif defined(LENCG)
				if(acc + rs2 + idxlens[j] < simT){
					accum[cid] = -2;
					#ifdef EXTRACOUNTS
					params->nPruneLength++;
					#endif
				}
				#endif
			} else if ( bnd >= simT && acc > -2 ){

	            if(qr == 0 && qdocs[cid] == 0){  // not a row we're interested in
	                acc = accum[cid] = -2;
                    cands[ncands++] = cid;
	            } else if ( acc == -1) {  // accumulate
					acc = accum[cid] = myval * idxvals[j];
					cands[ncands++] = cid;
				} else
					acc = accum[cid] += myval * idxvals[j];
				// check l2-norm bound
				#if defined(L2CG)
				if(acc + bsq * idxlens[j] < simTe){
					accum[cid] = -2;
					#ifdef EXTRACOUNTS
					params->nPruneLength++;
					#endif
				}
				#elif defined(LENCG)
				if(acc + rs2 + idxlens[j] < simT){
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
	da_stopwctimer(params->timer_1); // find candidates time

	da_startwctimer(params->timer_2); // process candidates time
#endif
	params->nCandidates  += ncands;
	params->nDotProducts += l2apProcessCandidates(params, rid, docs, lengths, sums, endptr, ncands,
			cands, accum, hashval, hashwgt, hashlen, hashsum, hashsz, rwgts, rfiwgts, ps);
#ifdef EXTRATIMES
	da_stopwctimer(params->timer_2); // process candidates time

	da_startwctimer(params->timer_1); // find candidates time
#endif
	// reset hashval
	for ( i=rowptr[rid]; i < endptr[rid]; i++ )
		hashval[rowind[i]] = 0;

#ifdef EXTRATIMES
	da_stopwctimer(params->timer_1); // find candidates time
#endif
}


/**
 * Reorder the document matrix in decreasing number of column nnz
 * and decreasing order of row max val. Also compute and store
 * necessary stats about the data
 * Store initial order in rperm/cperm.
 */
void l2apReorderDocs(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
		val_t *cwgts, idx_t *rperm, idx_t *cperm){
	ssize_t i, j, k, nnz;
	idx_t nrows, ncols;
	da_iikv_t *orderColNnz = NULL;
	da_ivkv_t *orderMaxVal = NULL;
	char pv = 0;
	ptr_t *rowptr, *nrowptr;
	idx_t *rowind, *nrowind;
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
	orderMaxVal = da_ivkvmalloc(nrows, "l2apReorderDocs: orderMaxVal");
	nrowptr = da_pmalloc(nrows+1, "l2apReorderDocs: nrowptr");
	nrowind = da_imalloc(nnz, "l2apReorderDocs: nrowind");
	nrowval = da_vmalloc(nnz, "l2apReorderDocs: nrowval");

	orderColNnz = da_iikvmalloc(ncols, "l2apReorderDocs: orderColNnz");
	//get new column order
	for(j=0; j < ncols; j++){
		orderColNnz[j].key = j;
		orderColNnz[j].val = csizes[j];
	}
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
	da_free((void**)&(*docs)->rowptr, &(*docs)->rowind, &(*docs)->rowval, LTERM);
	(*docs)->rowptr = nrowptr;
	(*docs)->rowind = nrowind;
	(*docs)->rowval = nrowval;

	da_csr_SortIndices((*docs), DA_ROW);

	//inverse the row permutation

	da_free((void**)&orderColNnz, &orderMaxVal, &cperm, LTERM);

	params->rperm = rperm;

}



/**
 * Reorder the document matrix in decreasing number of column nnzs.
 * Also compute and store necessary stats about the data.
 * Store initial order in cperm.
 */
void l2apReorderCols(da_csr_t *docs, idx_t *csizes, idx_t *cperm) {
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
    for (j=0; j < ncols; j++) {
        orderColNnz[j].key = j;
        orderColNnz[j].val = csizes[j];
    }
    da_iikvsortd(ncols, orderColNnz); // sort columns

    for (j=0; j < ncols; j++) {
        if(cperm)
            cperm[orderColNnz[j].key] = j;  // store column permutation
        csizes[j] = orderColNnz[j].val; // new column size after permuting
    }

    /* set new column ids in docs */
    for (i=0; i < nnz; i++)
        rowind[i] = cperm[rowind[i]];

    da_csr_SortIndices(docs, DA_ROW);

    da_free((void**) &orderColNnz, LTERM);

}
