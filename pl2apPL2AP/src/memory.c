/*!
 \file  memory.c
 \brief This file contains functions for allocating and freeing memory

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Aug 2013.

 \author David C. Anastasiu
 */

#include "includes.h"

/* This is for the global mcore that tracks all heap allocations */
static __thread da_mcore_t *damcore = NULL;


/*************************************************************************/
/*! This function allocates a two-dimensional matrix.
 */
/*************************************************************************/
void
da_AllocMatrix (void ***r_matrix, const size_t elmlen, const size_t ndim1,
		const size_t ndim2)
{
  size_t i, j;
  void **matrix;

  *r_matrix = NULL;

  if ((matrix =
       (void **) da_malloc (ndim1 * sizeof (void *),
			    "da_AllocMatrix: matrix")) == NULL)
    return;

  for (i = 0; i < ndim1; i++)
    {
      if ((matrix[i] =
	   (void *) da_malloc (ndim2 * elmlen,
			       "da_AllocMatrix: matrix[i]")) == NULL)
	{
	  for (j = 0; j < i; j++)
	    da_free ((void **) &matrix[j], LTERM);
	  return;
	}
    }

  *r_matrix = matrix;
}


/*************************************************************************/
/*! This function frees a two-dimensional matrix.
 */
/*************************************************************************/
void
da_FreeMatrix (void ***r_matrix, const size_t ndim1, const size_t ndim2)
{
  size_t i;
  void **matrix;

  if ((matrix = *r_matrix) == NULL)
    return;

  for (i = 0; i < ndim1; i++)
    da_free ((void **) &matrix[i], LTERM);

  da_free ((void **) r_matrix, LTERM);

}


/*************************************************************************/
/*! This function initializes tracking of heap allocations.
 */
/*************************************************************************/
int
da_malloc_init (void)
{
  if (damcore == NULL)
    damcore = da_gkmcoreCreate ();

  if (damcore == NULL)
    return 0;

  da_gkmcorePush (damcore);

  return 1;
}


/*************************************************************************/
/*! This function frees the memory that has been allocated since the
    last call to da_malloc_init().
 */
/*************************************************************************/
void
da_malloc_cleanup (const int showstats)
{
  if (damcore != NULL)
    {
      da_gkmcorePop (damcore);
      if (damcore->cmop == 0)
	{
	  da_gkmcoreDestroy (&damcore, showstats);
	  damcore = NULL;
	}
    }
}


/*************************************************************************/
/*! This function is my wrapper around malloc that provides the following
    enhancements over malloc:
 * It always allocates one byte of memory, even if 0 bytes are requested.
      This is to ensure that checks of returned values do not lead to NULL
      due to 0 bytes requested.
 * It zeros-out the memory that is allocated. This is for a quick init
      of the underlying datastructures.
 */
/**************************************************************************/
void *
da_malloc (const size_t nbytes, const char *const msg)
{
  void *ptr = NULL;

  if (nbytes > 0)
    {
      ptr = (void *) malloc (nbytes);
      if (ptr == NULL)
	{
	  fprintf (stderr, "   Current memory used:  %10zu bytes\n",
		   da_GetCurMemoryUsed ());
	  fprintf (stderr, "   Maximum memory used:  %10zu bytes\n",
		   da_GetMaxMemoryUsed ());
	  da_errexit
	    ("MEMORY ERROR: ***Memory allocation failed for %s. Requested size: %zu bytes",
	     msg, nbytes);
	  return NULL;
	}
      /* add this memory allocation */
      if (damcore != NULL)
	da_gkmcoreAdd (damcore, DA_MOPT_HEAP, nbytes, ptr);
    }
  else
    {
      ptr = (void *) malloc (1);
      if (ptr == NULL)
	{
	  fprintf (stderr, "   Current memory used:  %10zu bytes\n",
		   da_GetCurMemoryUsed ());
	  fprintf (stderr, "   Maximum memory used:  %10zu bytes\n",
		   da_GetMaxMemoryUsed ());
	  da_errexit
	    ("MEMORY ERROR: ***Memory allocation failed for %s. Requested size: 1 byte",
	     msg);
	  return NULL;
	}
      /* add this memory allocation */
      if (damcore != NULL)
	da_gkmcoreAdd (damcore, DA_MOPT_HEAP, 1, ptr);
    }

  return ptr;
}


/*************************************************************************
 * This function is my wrapper around realloc
 **************************************************************************/
void *
da_realloc (void *oldptr, const size_t nbytes, const char *const msg)
{
  void *ptr = NULL;

  /* remove this memory de-allocation */
  if (damcore != NULL && oldptr != NULL)
    da_gkmcoreDel (damcore, oldptr);

  if (nbytes > 0)
    {
      ptr = (void *) realloc (oldptr, nbytes);
      if (ptr == NULL)
	{
	  fprintf (stderr, "   Maximum memory used: %10zu bytes\n",
		   da_GetMaxMemoryUsed ());
	  fprintf (stderr, "   Current memory used: %10zu bytes\n",
		   da_GetCurMemoryUsed ());
	  da_errexit ("MEMORY ERROR: ***Memory realloc failed for %s. "
		      "Requested size: %zu bytes", msg, nbytes);
	  return NULL;
	}
      /* add this memory allocation */
      if (damcore != NULL)
	da_gkmcoreAdd (damcore, DA_MOPT_HEAP, nbytes, ptr);
    }
  else
    {
      ptr = (void *) realloc (oldptr, 1);
      if (ptr == NULL)
	{
	  fprintf (stderr, "   Maximum memory used: %10zu bytes\n",
		   da_GetMaxMemoryUsed ());
	  fprintf (stderr, "   Current memory used: %10zu bytes\n",
		   da_GetCurMemoryUsed ());
	  da_errexit ("MEMORY ERROR: ***Memory realloc failed for %s. "
		      "Requested size: %zu byte", msg);
	  return NULL;
	}
      /* add this memory allocation */
      if (damcore != NULL)
	da_gkmcoreAdd (damcore, DA_MOPT_HEAP, 1, ptr);
    }

  return ptr;
}


/*************************************************************************
 * This function is my wrapper around free, allows multiple pointers
 **************************************************************************/
void
da_free (void **ptr1, ...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
    {
      free (*ptr1);

      /* remove this memory de-allocation */
      if (damcore != NULL)
	da_gkmcoreDel (damcore, *ptr1);
    }
  *ptr1 = NULL;

  va_start (plist, ptr1);
  while ((ptr = va_arg (plist, void **)) != LTERM)
    {
      if (*ptr != NULL)
	{
	  free (*ptr);

	  /* remove this memory de-allocation */
	  if (damcore != NULL)
	    da_gkmcoreDel (damcore, *ptr);
	}
      *ptr = NULL;
    }
  va_end (plist);
}


/*************************************************************************
 * This function returns the current ammount of dynamically allocated
 * memory that is used by the system
 **************************************************************************/
size_t
da_GetCurMemoryUsed (void)
{
  if (damcore == NULL)
    return 0;
  else
    return damcore->cur_hallocs;
}


/*************************************************************************
 * This function returns the maximum ammount of dynamically allocated
 * memory that was used by the system
 **************************************************************************/
size_t
da_GetMaxMemoryUsed (void)
{
  if (damcore == NULL)
    return 0;
  else
    return damcore->max_hallocs;
}


/*************************************************************************/
/*! This function returns the VmSize and VmRSS of the calling process. */
/*************************************************************************/
void
da_GetVMInfo (size_t * vmsize, size_t * vmrss)
{
  FILE *fp;
  char fname[1024];

  sprintf (fname, "/proc/%d/statm", getpid ());
  fp = da_fopen (fname, "r", "proc/pid/statm");
  if (fscanf (fp, "%zu %zu", vmsize, vmrss) != 2)
    da_errexit ("Failed to read to values from %s\n", fname);
  da_fclose (fp);

  /*
   *vmsize *= sysconf(_SC_PAGESIZE);
   *vmrss  *= sysconf(_SC_PAGESIZE);
   */

  return;
}
