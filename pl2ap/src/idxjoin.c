/*!
 \file  idxjoin.c
 \brief This file contains idxJoin related functions

 The assumption in this implementation is that output fits in memory. Else we could not
 compute only against previous rows and expect output to be in sorted order.

 \author David C. Anastasiu
 */

#include "includes.h"


/**
 * Main entry point to IdxJoin.
 */
void
ijFindNeighbors (params_t * params)
{

  ssize_t i, j, k, nneighbs;
  idx_t nrows, ncands, progressInd;
  ptr_t *rowptr;
  idx_t *rowind, *marker = NULL;
  float eps;
  val_t *rowval;
  da_ivkv_t *hits = NULL, *cand = NULL;
  da_csr_t *docs, *neighbors = NULL;

  docs = params->docs;
  nrows = docs->nrows;		// num rows
  rowptr = docs->rowptr;	// index in rowind/rowval where each row starts
  rowind = docs->rowind;	// colid for each nnz
  rowval = docs->rowval;	// val for each nnz
  eps = params->eps;

  progressInd = nrows / 20;
  ncands = 0;			// number of considered candidates

  // allocate memory
  da_startwctimer (params->timer_5);	// memory allocation time

  simSearchSetup (params);

  hits = da_ivkvsmalloc (nrows, (da_ivkv_t)
			 {
			 0, 0.0}, "apFindNeighbors: hits");
  cand = da_ivkvsmalloc (nrows, (da_ivkv_t)
			 {
			 0, 0.0}, "apFindNeighbors: cand");
  marker = da_ismalloc (nrows, -1, "apFindNeighbors: marker");
  da_stopwctimer (params->timer_5);	// memory allocation time

  da_startwctimer (params->timer_7);	// indexing time
  da_csr_CreateIndex (docs, DA_COL);
  da_stopwctimer (params->timer_7);	// indexing time

  if (params->nim)
    {
#ifdef NO_OUTPUT
      da_errexit ("NO_OUTPUT directive invoked. Run without -nim.");
#endif

      neighbors = params->neighbors;

      // for each x \in V do
      da_startwctimer (params->timer_3);	// find neighbors time
      if (params->verbosity > 0)
	printf ("Progress Indicator (no. of records): ");
      for (i = 0; i < nrows; i++)
	{
	  k =
	    da_ApssGetSimilarSmallerRows (docs, i, 1,
					  rowptr[i + 1] - rowptr[i],
					  rowind + rowptr[i],
					  rowval + rowptr[i], DA_SIM_COS,
					  nrows, eps, hits, marker, cand,
					  &ncands);

	  params->nCandidates += ncands;
	  params->nDotProducts += ncands;

	  //transfer candidates
	  params->nSimPairs += k;
	  if (params->nim)
	    {
	      // set current num neighbors for this row's search
	      neighbors->rowptr[i + 1] = neighbors->rowptr[i];
	      // grow neighborhood matrix if necessary
	      if (neighbors->rowptr[i] + k > params->nghnnz)
		{
		  params->nghnnz += da_max (k, params->nghnnz / 2);
		  da_csr_Grow (neighbors, params->nghnnz);
		  params->nghinccnt++;
		}
	    }

	  for (j = 0; j < k; j++)
	    da_storeSim (params, i, hits[j].key, hits[j].val);

	  if (params->verbosity > 0 && i % progressInd == 0)
	    {
	      printf ("%d%%..", (int) ceil (100 * i / (float) nrows));
	      fflush (stdout);
	    }
	}
      if (params->verbosity > 0)
	printf ("100%%\n");
      da_stopwctimer (params->timer_3);	// find neighbors time

      simSearchFinalize (params);
    }
  else
    {

      // for each x \in V do
      da_startwctimer (params->timer_3);	// find neighbors time
      if (params->verbosity > 0)
	printf ("Progress Indicator (no. of records): ");
      for (i = 0; i < nrows; i++)
	{
//            if(params->symout)
	  k = da_ApssGetSimilarRows (docs, i, 1, rowptr[i + 1] - rowptr[i],
				     rowind + rowptr[i], rowval + rowptr[i],
				     DA_SIM_COS, nrows, eps, hits, marker,
				     cand, &ncands);
//            else
//                k = da_ApssGetSimilarSmallerRows(docs, i, 1, rowptr[i+1]-rowptr[i],
//                       rowind + rowptr[i], rowval + rowptr[i], DA_SIM_COS,
//                       nrows, eps, hits, marker, cand, &ncands);

	  params->nCandidates += ncands;
	  params->nDotProducts += ncands;

	  /* write the results in the file */
	  if (params->fpout)
	    {
	      if (params->verbosity > 4)
		printf (" writing... \n");
	      if (params->fmtWrite == DA_FMT_CSR)
		{
		  for (j = 0; j < k; j++)
		    {
		      fprintf (params->fpout, " %d %.7g", hits[j].key + 1,
			       hits[j].val);
		    }
		  fprintf (params->fpout, "\n");
		}
	      else
		for (j = 0; j < k; j++)
		  {
		    fprintf (params->fpout, "%zu %d %.7g\n", i + 1,
			     hits[j].key + 1, hits[j].val);
		  }
	      fflush (params->fpout);
	    }

	  params->nSimPairs += k;

	  if (params->verbosity > 0 && i % progressInd == 0)
	    {
	      printf ("%d%%..", (int) ceil (100 * i / (float) nrows));
	      fflush (stdout);
	    }
	}
      if (params->verbosity > 0)
	printf ("100%%\n");
      da_stopwctimer (params->timer_3);	// find neighbors time
//        if(params->symout)
      params->nSimPairs /= 2;
    }
  da_free ((void **) &hits, &cand, &marker, LTERM);
}



/**
 * Find similar rows in the matrix with smaller row id -  this version of the function reports
 * the number of candidates/dot products that were considered in the search.
 */
idx_t
da_ApssGetSimilarSmallerRows (da_csr_t * mat, idx_t rid, char noSelfSim,
			      idx_t nqterms, idx_t * qind, val_t * qval,
			      char simtype, idx_t nsim, float minsim,
			      da_ivkv_t * hits, idx_t * i_marker,
			      da_ivkv_t * i_cand, idx_t * ncands)
{
  ssize_t i, ii, j, k;
  idx_t nrows, ncols, ncand;
  ptr_t *colptr;
  idx_t *colind, *marker;
  val_t *colval, *rnorms, mynorm, *rsums, mysum;
  da_ivkv_t *cand;

  if (nqterms == 0)
    return 0;

  nrows = mat->nrows;
  ncols = mat->ncols;
  colptr = mat->colptr;
  colind = mat->colind;
  colval = mat->colval;

  marker =
    (i_marker ? i_marker :
     da_ismalloc (nrows, -1, "da_csr_GetSimilarSmallerRows: marker"));
  cand =
    (i_cand ? i_cand :
     da_ivkvmalloc (nrows, "da_csr_GetSimilarSmallerRows: cand"));

  switch (simtype)
    {
    case DA_SIM_COS:
      for (ncand = 0, ii = 0; ii < nqterms; ii++)
	{
	  i = qind[ii];
	  if (i < ncols)
	    {
	      for (j = colptr[i]; j < colptr[i + 1]; j++)
		{
		  k = colind[j];
		  if (k > rid || (noSelfSim && k == rid))
		    continue;
		  if (marker[k] == -1)
		    {
		      cand[ncand].key = k;
		      cand[ncand].val = 0;
		      marker[k] = ncand++;
		    }
		  cand[marker[k]].val += colval[j] * qval[ii];
		}
	    }
	}
      break;

    case DA_SIM_JAC:
      for (ncand = 0, ii = 0; ii < nqterms; ii++)
	{
	  i = qind[ii];
	  if (i < ncols)
	    {
	      for (j = colptr[i]; j < colptr[i + 1]; j++)
		{
		  k = colind[j];
		  if (k > rid || (noSelfSim && k == rid))
		    continue;
		  if (marker[k] == -1)
		    {
		      cand[ncand].key = k;
		      cand[ncand].val = 0;
		      marker[k] = ncand++;
		    }
		  cand[marker[k]].val += colval[j] * qval[ii];
		}
	    }
	}

      rnorms = mat->rnorms;
      mynorm = da_vdot (nqterms, qval, 1, qval, 1);

      for (i = 0; i < ncand; i++)
	cand[i].val =
	  cand[i].val / (rnorms[cand[i].key] + mynorm - cand[i].val);
      break;


    default:
      da_errexit ("Unknown similarity measure %d\n", simtype);
      return -1;
    }

  *ncands = ncand;

  /* go and prune the hits that are bellow minsim */
  for (j = 0, i = 0; i < ncand; i++)
    {
      marker[cand[i].key] = -1;
      if (cand[i].val >= minsim)
	cand[j++] = cand[i];
    }
  ncand = j;

  if (nsim == -1 || nsim >= ncand)
    {
      nsim = ncand;
    }
  else
    {
      nsim = da_min (nsim, ncand);
      da_ivkvkselectd (ncand, nsim, cand);
      da_ivkvsortd (nsim, cand);
    }

  da_ivkvcopy (nsim, cand, hits);

  if (i_marker == NULL)
    da_free ((void **) &marker, LTERM);
  if (i_cand == NULL)
    da_free ((void **) &cand, LTERM);

  return nsim;
}



/**
 * Find similar rows in the matrix -  this version of the function reports
 * the number of candidates/dot products that were considered in the search.
 */
idx_t
da_ApssGetSimilarRows (da_csr_t * mat, idx_t rid, char noSelfSim,
        idx_t nqterms, idx_t * qind, val_t * qval,
        char simtype, idx_t nsim, float minsim,
        da_ivkv_t * hits, idx_t * i_marker, da_ivkv_t * i_cand,
        idx_t * ncands)
{
    ssize_t i, ii, j, k;
    idx_t nrows, ncols, ncand;
    ptr_t *colptr;
    idx_t *colind, *marker;
    val_t *colval, *rnorms, mynorm, *rsums, mysum;
    da_ivkv_t *cand;

    if (nqterms == 0)
        return 0;

    nrows = mat->nrows;
    ncols = mat->ncols;
    colptr = mat->colptr;
    colind = mat->colind;
    colval = mat->colval;

    marker = (i_marker ? i_marker :
                    da_ismalloc (nrows, -1, "da_csr_GetSimilarSmallerRows: marker"));
    cand = (i_cand ? i_cand :
                    da_ivkvmalloc (nrows, "da_csr_GetSimilarSmallerRows: cand"));

    switch (simtype) {
        case DA_SIM_COS:
            for (ncand = 0, ii = 0; ii < nqterms; ii++) {
                i = qind[ii];
                if (i < ncols) {
                    for (j = colptr[i]; j < colptr[i + 1]; j++) {
                        k = colind[j];
                        if (noSelfSim && k == rid)
                            continue;
                        if (marker[k] == -1) {
                            cand[ncand].key = k;
                            cand[ncand].val = 0;
                            marker[k] = ncand++;
                        }
                        cand[marker[k]].val += colval[j] * qval[ii];
                    }
                }
            }
            break;

        case DA_SIM_JAC:
            for (ncand = 0, ii = 0; ii < nqterms; ii++) {
                i = qind[ii];
                if (i < ncols) {
                    for (j = colptr[i]; j < colptr[i + 1]; j++) {
                        k = colind[j];
                        if (noSelfSim && k == rid)
                            continue;
                        if (marker[k] == -1) {
                            cand[ncand].key = k;
                            cand[ncand].val = 0;
                            marker[k] = ncand++;
                        }
                        cand[marker[k]].val += colval[j] * qval[ii];
                    }
                }
            }

            rnorms = mat->rnorms;
            mynorm = da_vdot (nqterms, qval, 1, qval, 1);

            for (i = 0; i < ncand; i++)
                cand[i].val =
                        cand[i].val / (rnorms[cand[i].key] + mynorm - cand[i].val);
            break;


        default:
            da_errexit ("Unknown similarity measure %d\n", simtype);
            return -1;
    }

    *ncands = ncand;

    /* go and prune the hits that are bellow minsim */
    for (j = 0, i = 0; i < ncand; i++) {
        marker[cand[i].key] = -1;
        if (cand[i].val >= minsim)
            cand[j++] = cand[i];
    }
    ncand = j;

    if (nsim == -1 || nsim >= ncand) {
        nsim = ncand;
    } else {
        nsim = da_min (nsim, ncand);
        da_ivkvkselectd (ncand, nsim, cand);
        da_ivkvsortd (nsim, cand);
    }

    da_ivkvcopy (nsim, cand, hits);

    if (i_marker == NULL)
        da_free ((void **) &marker, LTERM);
    if (i_cand == NULL)
        da_free ((void **) &cand, LTERM);

    return nsim;
}
