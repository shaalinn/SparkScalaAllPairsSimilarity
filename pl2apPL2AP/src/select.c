/*!
\file  select.c
\brief Sorts only the largest k values

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Aug 2013.

\date   Started 2/28/2013
\author David C. Anastasiu
 */


#include "includes.h"

/* Byte-wise swap two items of size SIZE. */
#define DA_QSSWAP(a, b, stmp) do { stmp = (a); (a) = (b); (b) = stmp; } while (0)

/**** kselect functions for ix type kv arrays  ****/

/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t
da_ivkvkselectd (size_t n, idx_t topk, da_ivkv_t * cand)
{
  idx_t i, j, lo, hi, mid;
  da_ivkv_t stmp;
  val_t pivot;

  if (n <= topk)
    return n;			/* return if the array has fewer elements than we want */

  for (lo = 0, hi = n - 1; lo < hi;)
    {
      mid = lo + ((hi - lo) >> 1);

      /* select the median */
      if (cand[lo].val < cand[mid].val)
	mid = lo;
      if (cand[hi].val > cand[mid].val)
	mid = hi;
      else
	goto jump_over;
      if (cand[lo].val < cand[mid].val)
	mid = lo;

    jump_over:
      DA_QSSWAP (cand[mid], cand[hi], stmp);
      pivot = cand[hi].val;

      /* the partitioning algorithm */
      for (i = lo - 1, j = lo; j < hi; j++)
	{
	  if (cand[j].val >= pivot)
	    {
	      i++;
	      DA_QSSWAP (cand[i], cand[j], stmp);
	    }
	}
      i++;
      DA_QSSWAP (cand[i], cand[hi], stmp);


      if (i > topk)
	hi = i - 1;
      else if (i < topk)
	lo = i + 1;
      else
	break;
    }

  return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t
da_ivkvkselecti (size_t n, idx_t topk, da_ivkv_t * cand)
{
  idx_t i, j, lo, hi, mid;
  da_ivkv_t stmp;
  val_t pivot;

  if (n <= topk)
    return n;			/* return if the array has fewer elements than we want */

  for (lo = 0, hi = n - 1; lo < hi;)
    {
      mid = lo + ((hi - lo) >> 1);

      /* select the median */
      if (cand[lo].val > cand[mid].val)
	mid = lo;
      if (cand[hi].val < cand[mid].val)
	mid = hi;
      else
	goto jump_over;
      if (cand[lo].val > cand[mid].val)
	mid = lo;

    jump_over:
      DA_QSSWAP (cand[mid], cand[hi], stmp);
      pivot = cand[hi].val;

      /* the partitioning algorithm */
      for (i = lo - 1, j = lo; j < hi; j++)
	{
	  if (cand[j].val <= pivot)
	    {
	      i++;
	      DA_QSSWAP (cand[i], cand[j], stmp);
	    }
	}
      i++;
      DA_QSSWAP (cand[i], cand[hi], stmp);


      if (i > topk)
	hi = i - 1;
      else if (i < topk)
	lo = i + 1;
      else
	break;
    }

  return topk;
}



/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t
da_iakvkselectd (size_t n, idx_t topk, da_iakv_t * cand)
{
  idx_t i, j, lo, hi, mid;
  da_iakv_t stmp;
  val_t pivot;

  if (n <= topk)
    return n;			/* return if the array has fewer elements than we want */

  for (lo = 0, hi = n - 1; lo < hi;)
    {
      mid = lo + ((hi - lo) >> 1);

      /* select the median */
      if (cand[lo].val < cand[mid].val)
	mid = lo;
      if (cand[hi].val > cand[mid].val)
	mid = hi;
      else
	goto jump_over;
      if (cand[lo].val < cand[mid].val)
	mid = lo;

    jump_over:
      DA_QSSWAP (cand[mid], cand[hi], stmp);
      pivot = cand[hi].val;

      /* the partitioning algorithm */
      for (i = lo - 1, j = lo; j < hi; j++)
	{
	  if (cand[j].val >= pivot)
	    {
	      i++;
	      DA_QSSWAP (cand[i], cand[j], stmp);
	    }
	}
      i++;
      DA_QSSWAP (cand[i], cand[hi], stmp);


      if (i > topk)
	hi = i - 1;
      else if (i < topk)
	lo = i + 1;
      else
	break;
    }

  return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t
da_iakvkselecti (size_t n, idx_t topk, da_iakv_t * cand)
{
  idx_t i, j, lo, hi, mid;
  da_iakv_t stmp;
  val_t pivot;

  if (n <= topk)
    return n;			/* return if the array has fewer elements than we want */

  for (lo = 0, hi = n - 1; lo < hi;)
    {
      mid = lo + ((hi - lo) >> 1);

      /* select the median */
      if (cand[lo].val > cand[mid].val)
	mid = lo;
      if (cand[hi].val < cand[mid].val)
	mid = hi;
      else
	goto jump_over;
      if (cand[lo].val > cand[mid].val)
	mid = lo;

    jump_over:
      DA_QSSWAP (cand[mid], cand[hi], stmp);
      pivot = cand[hi].val;

      /* the partitioning algorithm */
      for (i = lo - 1, j = lo; j < hi; j++)
	{
	  if (cand[j].val <= pivot)
	    {
	      i++;
	      DA_QSSWAP (cand[i], cand[j], stmp);
	    }
	}
      i++;
      DA_QSSWAP (cand[i], cand[hi], stmp);


      if (i > topk)
	hi = i - 1;
      else if (i < topk)
	lo = i + 1;
      else
	break;
    }

  return topk;
}
