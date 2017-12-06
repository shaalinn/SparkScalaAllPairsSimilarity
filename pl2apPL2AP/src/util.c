/*!
 \file  util.c
 \brief This file contains utilities for the application

 \author David C. Anastasiu
 */

#include "includes.h"

/** forward declarations **/


/*************************************************************************/
/*! This function prints an error message and raises a signum signal
 */
/*************************************************************************/
void
da_errexit (const char *const f_str, ...)
{
    va_list argp;

    va_start (argp, f_str);
    vfprintf (stderr, f_str, argp);
    va_end (argp);

    fprintf (stderr, "\n");
    fflush (stderr);

    raise (SIGTERM);
}




/*************************************************************************
 * This function returns the CPU seconds
 **************************************************************************/
double
da_WClockSeconds (void)
{
#ifdef __GNUC__
    struct timeval ctime;

    gettimeofday (&ctime, NULL);

    return (double) ctime.tv_sec + (double) .000001 *ctime.tv_usec;
#else
    return (double) time (NULL);
#endif
}


/*************************************************************************
 * This function returns the CPU seconds
 **************************************************************************/
double
da_CPUSeconds (void)
{
    //#ifdef __OPENMP__
#ifdef __OPENMPXXXX__
    return omp_get_wtime ();
#else
#if defined(WIN32) || defined(__MINGW32__)
    return ((double) clock () / CLOCKS_PER_SEC);
#else
    struct rusage r;

    getrusage (RUSAGE_SELF, &r);
    return ((r.ru_utime.tv_sec + r.ru_stime.tv_sec) +
            1.0e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec));
#endif
#endif
}


/*************************************************************************/
/*! If format not specifically given (> 0), check if a text (non-binary) text file
 *  containing a csr is in CLUTO or CSR format.
    \param file is the matrix file to be checked.
    \return the CSR format: DA_FMT_CLUTO or DA_FMT_CSR
 */
/*************************************************************************/
char da_getFileFormat (char *file, const char format)
{
    if (format > 0)
        return format;
    size_t nnz;
    char fmt;
    char *ext, *p;

    ext = strrchr (file, '.');
    if (ext) {
        ext++;
        //make lowercase
        for (p = ext; *p; ++p)
            *p = tolower (*p);
        if ((fmt = da_getStringID (fmt_options, ext)) > -1)
            return fmt;
    } else if (da_fexists (file)) {
        // assume some sort of CSR. Can we guess?
        da_getfilestats (file, NULL, &nnz, NULL, NULL);
        return (nnz % 2 == 1) ? DA_FMT_CLUTO : DA_FMT_CSR;
    }
    return -1;
}

/**
 * Retrieve the dataset from the input filename
 */
char * da_getDataset (params_t * params)
{
    if (params->dataset)
        return params->dataset;
    int len;
    char *end, *start = strrchr (params->iFile, '/');
    if (start)
        start++;
    else
        start = params->iFile;
    end = strstr (start, ".");
    if (!end)
        end = params->iFile + strlen (params->iFile);
    len = end - start;
    params->dataset = da_cmalloc (len + 1, "da_getDataset: dataset");
    strncpy (params->dataset, start, len);
    params->dataset[len] = '\0';
    return params->dataset;
}


void da_gWriteVector (char *filename, val_t * vec, idx_t size, char *separator)
{
    idx_t i;
    FILE *fpout;

    if (filename)
        fpout = da_fopen (filename, "w", "da_vWriteVector: fpout");
    else
        fpout = stdout;

    for (i = 0; i < size; i++)
        fprintf (fpout, "%g%s", vec[i], separator);

    if (filename)
        da_fclose (fpout);
}


void da_WriteVector (char *filename, val_t * vec, idx_t size, char *separator)
{
    idx_t i;
    FILE *fpout;

    if (filename)
        fpout = da_fopen (filename, "w", "da_vWriteVector: fpout");
    else
        fpout = stdout;

    for (i = 0; i < size; i++)
        fprintf (fpout, PRNT_VALTYPE "%s", vec[i], separator);

    if (filename)
        da_fclose (fpout);
}

void da_vWriteMatrix (char *filename, val_t ** mat, idx_t nrows, idx_t ncols)
{
    idx_t i, j;
    FILE *fpout;

    if (filename)
        fpout = da_fopen (filename, "w", "da_vWriteMatrix: fpout");
    else
        fpout = stdout;

    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols - 1; j++)
            fprintf (fpout, PRNT_VALTYPE "\t", mat[i][j]);
        if (j < ncols)
            fprintf (fpout, PRNT_VALTYPE, mat[i][j]);
        fprintf (fpout, "\n");
    }

    if (filename)
        da_fclose (fpout);
}

void da_dWriteMatrix (char *filename, double **mat, idx_t nrows, idx_t ncols)
{
    idx_t i, j;
    FILE *fpout;

    if (filename)
        fpout = da_fopen (filename, "w", "da_vWriteMatrix: fpout");
    else
        fpout = stdout;

    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols - 1; j++)
            fprintf (fpout, "%f\t", mat[i][j]);
        if (j < ncols)
            fprintf (fpout, "%f", mat[i][j]);
        fprintf (fpout, "\n");
    }

    if (filename)
        da_fclose (fpout);
}

void da_iWriteMatrix (char *filename, idx_t ** mat, idx_t nrows, idx_t ncols)
{
    idx_t i, j;
    FILE *fpout;

    if (filename)
        fpout = da_fopen (filename, "w", "da_iWriteMatrix: fpout");
    else
        fpout = stdout;

    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols - 1; j++)
            fprintf (fpout, PRNT_IDXTYPE "\t", mat[i][j]);
        if (j < ncols)
            fprintf (fpout, PRNT_IDXTYPE, mat[i][j]);
        fprintf (fpout, "\n");
    }

    if (filename)
        da_fclose (fpout);
}

void da_iWriteVector (char *filename, idx_t * vec, idx_t size, char *separator)
{
    idx_t i;
    FILE *fpout;

    if (filename)
        fpout = da_fopen (filename, "w", "da_iWriteVector: fpout");
    else
        fpout = stdout;

    for (i = 0; i < size; i++)
        fprintf (fpout, PRNT_IDXTYPE "%s", vec[i], separator);

    if (filename)
        da_fclose (fpout);
}


void da_pWriteVector (char *filename, idx_t * vec, idx_t size, char *separator)
{
    idx_t i;
    FILE *fpout;

    if (filename)
        fpout = da_fopen (filename, "w", "da_iWriteVector: fpout");
    else
        fpout = stdout;

    for (i = 0; i < size; i++)
        fprintf (fpout, PRNT_IDXTYPE "%s", vec[i], separator);

    if (filename)
        da_fclose (fpout);
}

/*************************************************************************/
/*! Stores similarities in the neighbors matrix or writes them out to a file
 * in ijv format. Uses cache to buffer writes. Empties buffer when i is negative.
 *
 *	\param fpout file to write output to
 *	\param cache structure to hold similarity values temporarily
 *	\param zidx Starting col index identifier (1 if first col is 1)
 *	\param i the query neighbor
 *	\param j the other of the two neighbors
 *	\param v similarity value
 */
/*************************************************************************/
void da_storeSim (params_t * params, idx_t i, idx_t j, val_t v)
{

#ifndef NO_OUTPUT
    FILE *fpout;
    da_sims_t *cache;
    idx_t zidx;
    int32_t k;
    da_sim_t *sims;
    ssize_t nneighbs;

    if (params->nim) {
        /* we write neighbors in memory for sake of experiments, even if not saving to stdout or disk */
        nneighbs = params->neighbors->rowptr[i + 1];
        params->neighbors->rowind[nneighbs] = j;
        params->neighbors->rowval[nneighbs] = v;
        params->neighbors->rowptr[i + 1]++;
    } else if (params->fpout) {

        fpout = params->fpout;
        cache = params->neighcache;
        zidx = params->writeNum;

        if (i >= 0) {
            if (params->rperm) {
                i = params->rperm[i];
                j = params->rperm[j];
            }
            cache->sims[cache->nsims].i = i;
            cache->sims[cache->nsims].j = j;
            cache->sims[cache->nsims++].val = v;
        }

        if (i < 0 || cache->nsims == cache->size) {
            sims = cache->sims;
            for (k = 0; k < cache->nsims; k++) {
                if (params->symout)
                    fprintf (fpout, "%d %d %.8f\n%d %d %.8f\n",
                            (int32_t) sims[k].i + zidx,
                            (int32_t) sims[k].j + zidx, (float) sims[k].val,
                            (int32_t) sims[k].j + zidx,
                            (int32_t) sims[k].i + zidx, (float) sims[k].val);
                else if (sims[k].i < sims[k].j)
                    fprintf (fpout, "%d %d %.8f\n", (int32_t) sims[k].i + zidx,
                            (int32_t) sims[k].j + zidx, (float) sims[k].val);
                else
                    fprintf (fpout, "%d %d %.8f\n", (int32_t) sims[k].j + zidx,
                            (int32_t) sims[k].i + zidx, (float) sims[k].val);

            }
            cache->nsims = 0;
        }

    }
#endif
}

char da_isFmtBinary (char fmt)
{
    switch (fmt)
    {
        case DA_FMT_BIJV:
        case DA_FMT_BINAP:
        case DA_FMT_BINAPB:
        case DA_FMT_BINCOL:
        case DA_FMT_BINROW:
            return 1;
        default:
            return 0;
    }
}


/*************************************************************************/
/*! Prints a timer in human readable format. Timers should be started and
 * stopped with da_startwctimer(t) and da_stopwctimer(t) respectively.
 *
 *	\param name is the timer name
 *	\param time is the time in seconds
 */
/*************************************************************************/
void da_printTimer (char *name, double time)
{
    if (time < 60) {
        printf ("%s %.2fs\n", name, time);
    } else if (time < 3600) {
        int32_t min = time / 60;
        time -= min * 60;
        printf ("%s %dm %.2fs\n", name, min, time);
    } else {
        int32_t hours = time / 3600;
        time -= hours * 3600;
        int32_t min = time / 60;
        time -= min * 60;
        printf ("%s %dh %dm %.2fs\n", name, hours, min, time);
    }
}

/************************************************************************/
/*! \brief Prints a timer in long format, including human readable and
 * computer readable time. Timers should be started and stopped with
 * da_startwctimer(t) and da_stopwctimer(t) respectively.

 	\param name is the Timer name, a string to be printed before the timer.
 	\param time is the time that should be reported.
 */
/*************************************************************************/
void da_printTimerLong (char *name, double time)
{
    double t = time;
    if (time < 60) {
        printf ("%s %.2f (%.2fs)\n", name, t, time);
    }
    else if (time < 3600) {
        int32_t min = time / 60;
        time -= min * 60;
        printf ("%s %.2f (%dm %.2fs)\n", name, t, min, time);
    } else {
        int32_t hours = time / 3600;
        time -= hours * 3600;
        int32_t min = time / 60;
        time -= min * 60;
        printf ("%s %.2f (%dh %dm %.2fs)\n", name, t, hours, min, time);
    }
}


/**
 * http://www.codemaestro.com/reviews/9
 */
float da_fsqrt (float number)
{
    long i;
    float x, y;
    const float f = 1.5F;

    x = number * 0.5F;
    y = number;
    i = *(long *) &y;
    i = 0x5f3759df - (i >> 1);
    y = *(float *) &i;
    y = y * (f - (x * y * y));
    y = y * (f - (x * y * y));
    return number * y;
}

float da_finvSqrt (float x)
{
    float xhalf = 0.5f * x;
    int i = *(int *) &x;		// get bits for floating value
    i = 0x5f375a86 - (i >> 1);	// gives initial guess y0
    x = *(float *) &i;		// convert bits back to float
    x = x * (1.5f - xhalf * x * x);	// Newton step, repeating increases accuracy
    return x;
}



/**
 * Add the transpose of a matrix to its row structure, ignoring self-similarity
 * NOTE1: Assumes rowind has sorted columns
 * NOTE2: Method allocates nnz*2 memory for nval & nrowval.
 *        New matrix may have less nnzs. One may wish to check and
 *        realloc these arrays to free unused space. We leave this
 *        as an optional external procedure.
 */
void da_addHigherNeighbors (da_csr_t * mat)
{
    ssize_t i, j, k, l, nnz, nrows;
    ptr_t *nrowptr = NULL, *rowptr = NULL, *colptr = NULL;
    idx_t *rowind = NULL, *colind = NULL, *nrowind = NULL, *itmp = NULL;
    val_t *nrowval = NULL, *rowval = NULL, *colval = NULL, *ftmp = NULL;

    ASSERT (mat->nrows == mat->ncols);

    if (mat->nrows < 1)
        return;

    nnz = mat->rowptr[mat->nrows];
    nrows = mat->nrows;
    rowptr = mat->rowptr;
    rowind = mat->rowind;
    rowval = mat->rowval;
    nrowptr = da_pmalloc (nrows + 1, "da_addHigherNeighbors: nrowptr");
    nrowind = da_imalloc (nnz * 2, "da_addHigherNeighbors: nrowind");
    nrowval = da_vmalloc (nnz * 2, "da_addHigherNeighbors: nrowval");
    nrowptr[0] = 0;

    // create the transpose of the matrix in the col strucure of mat
    da_csr_CreateIndex (mat, DA_COL);
    da_csr_SortIndices (mat, DA_COL);
    // by definition, the col structure of the mat is equivalent to the transpose of the
    // matrix if thinking of it as row structure instead of col structure
    colptr = mat->colptr;
    colind = mat->colind;
    colval = mat->colval;

    // add matrices in row and col structures & store in nrow structure
    for (i = 0, nnz = 0; i < nrows; i++) {
        for (j = rowptr[i], k = colptr[i];
                j < rowptr[i + 1] && k < colptr[i + 1];) {
            if (rowind[j] == colind[k]) {
                nrowval[nnz] = da_max (rowval[j], colval[k]);
                nrowind[nnz++] = rowind[j];
                j++;
                k++;
            } else if (rowind[j] < colind[k]) {
                nrowval[nnz] = rowval[j];
                nrowind[nnz++] = rowind[j];
                j++;
            } else {
                nrowval[nnz] = colval[k];
                nrowind[nnz++] = colind[k];
                k++;
            }
        }
        for (; j < rowptr[i + 1]; j++) {
            nrowval[nnz] = rowval[j];
            nrowind[nnz++] = rowind[j];
        }
        for (; k < colptr[i + 1]; k++) {
            nrowval[nnz] = colval[k];
            nrowind[nnz++] = colind[k];
        }
        nrowptr[i + 1] = nnz;
    }

    da_free ((void **) &mat->rowptr, &mat->rowind, &mat->rowval,
            &mat->colptr, &mat->colind, &mat->colval, LTERM);
    mat->rowptr = nrowptr;
    mat->rowind = nrowind;
    mat->rowval = nrowval;

}


/*************************************************************************
 * Combine the output matrices from several search results and undo the permutation
 * Perm contains the original ids associated with rows in docs matrix
 **************************************************************************/
da_csr_t * da_csr_CombineAndPermute (idx_t nrows, idx_t sz, da_tsearch_t ** tsearch,
        idx_t * perm)
{
    size_t i, j, k, ii, jj, nnz, orgI, orgJ;
    da_tsearch_t *ts;
    da_csr_t *nbrs, *mat = NULL;
    idx_t *ind;
    ptr_t *ptr;
    val_t *val;

    if (sz < 1)
        return NULL;
    ptr = da_psmalloc (nrows + 1, 0, "da_csr_CombineAndPermute: cnt");

    mat = da_csr_Create ();
    mat->nrows = mat->ncols = nrows;
    mat->rowptr = ptr =
            da_psmalloc (nrows + 1, 0, "da_csr_CombineAndPermute: mat->rowperm");
    for (nnz = 0, k = 0; k < sz; ++k) {
        ts = tsearch[k];
        nbrs = ts->nbrs;
        nnz += nbrs->rowptr[nbrs->nrows];
        for (i = 0; i < nbrs->nrows; ++i) {
            ii = ts->nbrids[i];	/** id in docs for query row **/
            orgI = perm[ii];	/** id in original order for query row **/
            ptr[orgI] = nbrs->rowptr[i + 1] - nbrs->rowptr[i];
        }
    }
    CSRMAKE (i, nrows, ptr);
    mat->rowind = ind =
            da_ismalloc (nnz, 0, "da_csr_CombineAndPermute: mat->rowind");
    mat->rowval = val =
            da_vsmalloc (nnz, 0, "da_csr_CombineAndPermute: mat->rowval");


    for (nnz = 0, k = 0; k < sz; ++k) {
        ts = tsearch[k];
        nbrs = ts->nbrs;
        for (i = 0; i < nbrs->nrows; ++i) {
            ii = ts->nbrids[i];
            orgI = perm[ii];
            for (jj = ptr[orgI], j = nbrs->rowptr[i]; j < nbrs->rowptr[i + 1];
                    ++j, ++jj) {
                ind[jj] = perm[nbrs->rowind[j]];
                val[jj] = nbrs->rowval[j];
            }
        }

    }

    return mat;
}


/*************************************************************************
 * Combine the output matrices from several search results and undo the permutation
 * Perm contains the original ids associated with rows in docs matrix
 * Items are re-arranged to fit in lower-triangular section of the output matrix
 **************************************************************************/
da_csr_t * da_csr_CombineAndPermuteLT (idx_t nrows, idx_t sz, da_tsearch_t ** tsearch,
        idx_t * perm)
{
    size_t i, j, k, ii, jj, nnz, orgI, orgII, orgJ;
    da_tsearch_t *ts;
    da_csr_t *nbrs, *mat = NULL;
    da_sim_t *nbrslst;
    idx_t *ind;
    ptr_t *ptr;
    val_t *val;

    if (sz < 1)
        return NULL;

    nnz = 0;
    mat = da_csr_Create ();
    mat->nrows = mat->ncols = nrows;
    mat->rowptr = ptr =
            da_psmalloc (nrows + 1, 0, "da_csr_CombineAndPermute: mat->rowptr");
    if (tsearch[0]->nbrs) {
        for (nnz = 0, k = 0; k < sz; ++k) {
            ts = tsearch[k];
            nbrs = ts->nbrs;
            nnz += nbrs->rowptr[nbrs->nrows];
            for (i = 0; i < nbrs->nrows; ++i) {
                ii = ts->nbrids[i];   /** id in docs for query row **/
                orgI = perm[ii];	    /** id in original order for query row **/
                for (j = nbrs->rowptr[i]; j < nbrs->rowptr[i + 1]; ++j) {
                    orgII = perm[nbrs->rowind[j]];
                    if (orgI > orgII) {
                        ptr[orgI]++;
                    } else {
                        ptr[orgII]++;
                    }
                }
            }
        }
    }
    else if (tsearch[0]->nbrslst) {
        for (nnz = 0, k = 0; k < sz; ++k) {
            ts = tsearch[k];
            nbrslst = ts->nbrslst;
            nnz += ts->nbrnnz;
            for (i = 0; i < ts->nbrnnz; ++i) {
                orgI = perm[nbrslst[i].i];      /** id in original order for query row **/
                orgII = perm[nbrslst[i].j];     /** id in original order for candidate row **/
                if (orgI > orgII) {
                    ptr[orgI]++;
                } else {
                    ptr[orgII]++;
                }
            }
        }
    }
    else
        da_errexit ("Could not find any neighbors.");
    CSRMAKE (i, nrows, ptr);
    mat->rowind = ind =
            da_ismalloc (nnz, 0, "da_csr_CombineAndPermute: mat->rowind");
    mat->rowval = val =
            da_vsmalloc (nnz, 0, "da_csr_CombineAndPermute: mat->rowval");

    if (tsearch[0]->nbrs) {
        for (nnz = 0, k = 0; k < sz; ++k) {
            ts = tsearch[k];
            nbrs = ts->nbrs;
            for (i = 0; i < nbrs->nrows; ++i) {
                ii = ts->nbrids[i];
                orgI = perm[ii];
                for (j = nbrs->rowptr[i]; j < nbrs->rowptr[i + 1]; ++j) {
                    orgII = perm[nbrs->rowind[j]];
                    if (orgI > orgII) {
                        jj = ptr[orgI]++;
                        ind[jj] = orgII;
                        val[jj] = nbrs->rowval[j];
                    } else {
                        jj = ptr[orgII]++;
                        ind[jj] = orgI;
                        val[jj] = nbrs->rowval[j];
                    }
                }
            }
        }
    } else if (tsearch[0]->nbrslst) {
        for (nnz = 0, k = 0; k < sz; ++k) {
            ts = tsearch[k];
            nbrslst = ts->nbrslst;
            for (i = 0; i < ts->nbrnnz; ++i) {
                orgI = perm[nbrslst[i].i];      /** id in original order for query row **/
                orgII = perm[nbrslst[i].j];     /** id in original order for candidate row **/
                if (orgI > orgII) {
                    jj = ptr[orgI]++;
                    ind[jj] = orgII;
                    val[jj] = nbrslst[i].val;
                } else {
                    jj = ptr[orgII]++;
                    ind[jj] = orgI;
                    val[jj] = nbrslst[i].val;
                }
            }
        }
    }
    CSRSHIFT (i, nrows, ptr);

    return mat;
}



/**
 * Undo the permutation of the output square matrix, storing results in lower-triangular form
 */
void da_inversePermuteSqMatrixLT (da_csr_t * mat, idx_t * perm, char which)
{
    ssize_t i, j, jj, orgI, orgII, orgJ, nrows;
    ptr_t nnz;
    ptr_t *ptr = NULL, *nptr = NULL;
    idx_t *ind = NULL, *nind = NULL;
    val_t *val = NULL, *nval = NULL;

    if (!mat->rowptr)
        da_errexit
        ("da_inversePermuteSqMatrixLT: Row structure of the matrix is missing.");

    nrows = mat->nrows;
    nnz = mat->rowptr[nrows];
    ptr = mat->rowptr;
    ind = mat->rowind;
    val = mat->rowval;
    nptr = da_psmalloc (nrows + 1, 0, "da_inversePermuteSqMatrixLT: nptr");
    nind = da_ismalloc (nnz, 0, "da_inversePermuteSqMatrixLT: nind");
    nval = da_vsmalloc (nnz, 0.0, "da_inversePermuteSqMatrixLT: nval");

    for (i = 0; i < nrows; i++)
    {
        orgI = perm[i];
        for (j = ptr[i]; j < ptr[i + 1]; ++j) {
            orgII = perm[ind[j]];
            ASSERT (orgI < nrows && orgII < nrows);
            if (orgI > orgII) {
                nptr[orgI]++;
            } else {
                nptr[orgII]++;
            }
        }
    }
    CSRMAKE (i, nrows, nptr);

    for (i = 0; i < nrows; i++) {
        orgI = perm[i];
        for (j = ptr[i]; j < ptr[i + 1]; ++j) {
            orgII = perm[ind[j]];
            ASSERT (orgI < nrows && orgII < nrows);
            if (orgI > orgII) {
                jj = nptr[orgI]++;
                nind[jj] = (idx_t) orgII;
                nval[jj] = val[j];
            } else {
                jj = nptr[orgII]++;
                nind[jj] = (idx_t) orgI;
                nval[jj] = val[j];
            }
        }
    }
    CSRSHIFT (i, nrows, nptr);

    da_pcopy (nrows + 1, nptr, ptr);
    da_icopy (nnz, nind, ind);
    da_vcopy (nnz, nval, val);
    da_free ((void **) &nptr, &nind, &nval, LTERM);

    if (which & DA_ROW)
        da_csr_SortIndices (mat, DA_ROW);
    if (which & DA_COL)
        da_csr_SortIndices (mat, DA_COL);

}


/**
 * Undo the matrix permutation
 */
void da_inversePermuteMatrix (da_csr_t ** matP, idx_t * rowPerm, idx_t * colPerm)
{
    ssize_t i, j, k, orgI, orgJ, sz;
    ptr_t nnz;
    ptr_t *ptr, *nptr;
    idx_t *ind, *nind, *rperm, *cperm, *revperm = NULL, *itmp = NULL;
    val_t *val, *nval = NULL, *vtmp = NULL;
    da_csr_t *mat = *matP;

    if (mat->rowptr) {

        sz = mat->nrows;
        nnz = mat->rowptr[sz];
        ptr = mat->rowptr;
        ind = mat->rowind;
        val = mat->rowval;
        rperm = rowPerm;
        cperm = colPerm;
        nptr = da_pmalloc (sz + 1, "da_inversePermuteMatrix: nptr");
        nind = da_imalloc (nnz, "da_inversePermuteMatrix: nind");
        nval = da_vmalloc (nnz, "da_inversePermuteMatrix: nval");
        revperm = da_imalloc (sz, "da_inversePermuteMatrix: revperm");
        nptr[0] = 0;

        for (i = 0; i < sz; i++)
            revperm[rperm[i]] = i;

        for (i = 0; i < sz; i++) {
            orgI = revperm[i];
            nptr[i + 1] = nptr[i] + (ptr[orgI + 1] - ptr[orgI]);
            for (j = nptr[i]; j < nptr[i + 1]; j++) {
                orgJ = ptr[orgI] + j - nptr[i];
                nind[j] = cperm[ind[orgJ]];
                if (val)
                    nval[j] = val[orgJ];
            }
        }

        da_free ((void **) &(*matP)->rowptr, &(*matP)->rowind, &(*matP)->rowval,
                LTERM);
        (*matP)->rowptr = nptr;
        (*matP)->rowind = nind;
        (*matP)->rowval = nval;

        da_csr_SortIndices (mat, DA_ROW);

        da_free ((void **) &revperm, LTERM);
    }

    if (mat->colptr) {

        sz = mat->ncols;
        nnz = mat->colptr[sz];
        ptr = mat->colptr;
        ind = mat->colind;
        val = mat->colval;
        rperm = colPerm;
        cperm = rowPerm;
        nptr = da_pmalloc (sz + 1, "da_inversePermuteMatrix: nptr");
        nind = da_imalloc (nnz, "da_inversePermuteMatrix: nind");
        nval = da_vmalloc (nnz, "da_inversePermuteMatrix: nval");
        revperm = da_imalloc (sz, "da_inversePermuteMatrix: revperm");
        nptr[0] = 0;

        for (i = 0; i < sz; i++)
            revperm[rperm[i]] = i;

        for (i = 0; i < sz; i++) {
            orgI = revperm[i];
            nptr[i + 1] = nptr[i] + (ptr[orgI + 1] - ptr[orgI]);
            for (j = nptr[i]; j < nptr[i + 1]; j++) {
                orgJ = ptr[orgI] + j - nptr[i];
                nind[j] = cperm[ind[orgJ]];
                if (val)
                    nval[j] = val[orgJ];
            }
        }

        da_free ((void **) &(*matP)->colptr, &(*matP)->colind, &(*matP)->colval,
                LTERM);
        (*matP)->colptr = nptr;
        (*matP)->colind = nind;
        (*matP)->colval = nval;

        da_csr_SortIndices ((*matP), DA_COL);

        da_free ((void **) &revperm, LTERM);
    }


}



/**
 * Undo the matrix permutation in place
 */
void da_inversePermuteRows (da_csr_t * mat, idx_t * rowPerm, idx_t * colPerm)
{
    ssize_t i, j, k, orgI, orgJ, sz;
    ptr_t nnz;
    ptr_t *ptr, *nptr;
    idx_t *ind, *nind, *rperm, *cperm, *revperm = NULL, *itmp = NULL;
    val_t *val, *nval = NULL, *vtmp = NULL;

    if (!rowPerm)
        return;

    if (!mat->rowptr)
        da_errexit ("da_inversePermuteRows: row structure required!");

    sz = mat->nrows;
    nnz = mat->rowptr[sz];
    ptr = mat->rowptr;
    ind = mat->rowind;
    val = mat->rowval;
    rperm = rowPerm;
    cperm = colPerm;
    nptr = da_pmalloc (sz + 1, "da_inversePermuteMatrix: nptr");
    nind = da_imalloc (nnz, "da_inversePermuteMatrix: nind");
    nval = da_vmalloc (nnz, "da_inversePermuteMatrix: nval");
    revperm = da_imalloc (sz, "da_inversePermuteMatrix: revperm");
    nptr[0] = 0;

    for (i = 0; i < sz; i++)
        revperm[rperm[i]] = i;

    for (i = 0; i < sz; i++) {
        orgI = revperm[i];
        nptr[i + 1] = nptr[i] + (ptr[orgI + 1] - ptr[orgI]);
        for (j = nptr[i]; j < nptr[i + 1]; j++) {
            orgJ = ptr[orgI] + j - nptr[i];
            nind[j] = cperm[ind[orgJ]];
            if (val)
                nval[j] = val[orgJ];
        }
    }

    da_pcopy (sz + 1, nptr, ptr);
    da_icopy (nnz, nind, ind);
    if (val)
        da_vcopy (nnz, nval, val);

    da_csr_SortIndices (mat, DA_ROW);

    da_free ((void **) &revperm, &nptr, &nind, &nval, LTERM);

}



double da_CosineSimilarity (const size_t nind1, const idx_t const *ind1,
        const val_t const *val1, const size_t nind2,
        const idx_t const *ind2, const val_t const *val2)
{
    int i1, i2;
    double norm1, norm2;
    double sim;

    sim = 0.0;
    i1 = i2 = 0;
    while (i1 < nind1 && i2 < nind2) {
        if (ind1[i1] < ind2[i2]) {
            i1++;
        } else if (ind1[i1] > ind2[i2]) {
            i2++;
        } else {
            sim += val1[i1] * val2[i2];
            i1++;
            i2++;
        }
    }

    return sim;
}





/*************************************************************************
 * This function returns the log2(x)
 **************************************************************************/
int da_log2 (idx_t a)
{
    ssize_t i;

    for (i = 1; a > 1; i++, a = a >> 1);
    return i - 1;
}


/*************************************************************************
 * This function checks if the argument is a power of 2
 **************************************************************************/
int da_ispow2 (idx_t a)
{
    return (a == (1 << da_log2 (a)));
}


/*************************************************************************
 * This function returns the log2(x)
 **************************************************************************/
float da_flog2 (float a)
{
    return log (a) / log (2.0);
}

/*************************************************************************
 * This function returns the log2(x)
 **************************************************************************/
double da_dlog2 (double a)
{
    return log (a) / log (2.0);
}

/*************************************************************************
 * This function returns the log2(x)
 **************************************************************************/
val_t da_vlog2 (val_t a)
{
    return log (a) / log (2.0);
}




#define COMPERRPRINT(ar, ac, av, br, bc, bv, rc, fr) \
        printf("%sa[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_VALTYPE \
                    "] != b[" PRNT_IDXTYPE "," PRNT_IDXTYPE "," PRNT_VALTYPE "]", fr++?", ":"", rc?(ar):(ac), rc?(ac):(ar), av, rc?(br):(bc), rc?(bc):(br), bv);

#define COMPERRPRINT_A(ar, ac, bv, rc, fr) \
        printf("%s!a[" PRNT_IDXTYPE "," PRNT_IDXTYPE ",(" PRNT_VALTYPE ")]", fr++?", ":"", rc?(ar):(ac), rc?(ac):(ar), bv);

#define COMPERRPRINT_B(br, bc, av, rc, fr) \
        printf("%s!b[" PRNT_IDXTYPE "," PRNT_IDXTYPE ",(" PRNT_VALTYPE ")]", fr++?", ":"", rc?(br):(bc), rc?(bc):(br), av);

#define COMPERRPRINT_VAL(ar, ac, av, bv, rc, fr) \
        printf("%s![" PRNT_IDXTYPE "," PRNT_IDXTYPE ",(" PRNT_VALTYPE ", " PRNT_VALTYPE ")]", fr++?", ":"", rc?(ar):(ac), rc?(ac):(ar), av, bv);


/**
 * Compare two csr matrices and print out differences
 * 	\param doc1 first matrix to compare
 * 	\param doc2 second matrix to compare
 * 	\param eps Float max delta for value comparison
 * 	\param compVals Whether values should be compared
 */
void da_csrCompare (da_csr_t * doc1, da_csr_t * doc2, float eps, char compInds,
        char compVals)
{
    ssize_t j, k, ndiff = 0;
    idx_t i, l, fr;
    da_csr_t *a = NULL, *b = NULL;
    ptr_t *ptr1, *ptr2;
    idx_t *ind1, *ind2;
    val_t *val1, *val2;
    char rc;

    ASSERT ((doc1->rowptr && doc2->rowptr) || (doc1->colptr && doc2->colptr));

    a = da_csr_Copy (doc1);
    b = da_csr_Copy (doc2);
    if (compInds) {
        da_csr_SortIndices (a, DA_ROW);
        da_csr_SortIndices (b, DA_ROW);
    }

    if (a->rowptr) {
        rc = 1;
        ptr1 = a->rowptr;
        ind1 = a->rowind;
        val1 = a->rowval;
        ptr2 = b->rowptr;
        ind2 = b->rowind;
        val2 = b->rowval;
    } else {
        rc = 0;
        ptr1 = a->rowptr;
        ind1 = a->rowind;
        val1 = a->rowval;
        ptr2 = b->rowptr;
        ind2 = b->rowind;
        val2 = b->rowval;
    }
    if (!val1)
        compVals = 0;

    if ((a->nrows != b->nrows) || (a->ncols != b->ncols)
            || (ptr1[a->nrows] != ptr2[b->nrows]))
        printf ("Matrix stats differ: A[" PRNT_IDXTYPE "," PRNT_IDXTYPE ","
                PRNT_PTRTYPE "] != B[" PRNT_IDXTYPE "," PRNT_IDXTYPE ","
                PRNT_PTRTYPE "].\n", a->nrows, a->ncols, ptr1[a->nrows], b->nrows,
                b->ncols, ptr2[b->nrows]);
    printf ("Differences: \n");
    if (compVals && !compInds) {
        for (i = 0; i < a->nrows && i < b->nrows; i++) {
            fr = 0;
            for (j = ptr1[i], k = ptr2[i]; j < ptr1[i + 1] && k < ptr2[i + 1];
                    ++j, ++k) {
                if (da_abs (val1[j] - val2[k]) > eps) {
                    printf ("%s[%d %zu %f %f]", fr++ ? ", " : "", i,
                            j - ptr1[i] + 1, val1[j], val2[k]);
                    ndiff++;
                }
            }
            for (; j < ptr1[i + 1]; j++) {
                printf ("%s!b[%d %zu %f]", fr++ ? ", " : "", i, j - ptr1[i] + 1,
                        val1[j]);
                ndiff++;
            }
            for (; k < ptr2[i + 1]; k++) {
                printf ("%s!a[%d %zu %f]", fr++ ? ", " : "", i, k - ptr2[i] + 1,
                        val2[k]);
                ndiff++;
            }
            if (fr)
                printf ("\n");
        }
        for (l = i; i < a->nrows; i++) {
            for (fr = 0, j = ptr1[i]; j < ptr1[i + 1]; j++) {
                printf ("%s!b[%d %zu %f]", fr++ ? ", " : "", i, j - ptr1[i] + 1,
                        val1[j]);
                ndiff++;
            }
            if (fr)
                printf ("\n");
        }
        i = l;
        for (; i < b->nrows; i++) {
            for (fr = 0, k = ptr2[i]; k < ptr2[i + 1]; k++) {
                printf ("%s!a[%d %zu %f]", fr++ ? ", " : "", i, k - ptr2[i] + 1,
                        val2[k]);
                ndiff++;
            }
            if (fr)
                printf ("\n");
        }

    } else {

        for (i = 0; i < a->nrows && i < b->nrows; i++) {
            fr = 0;
            for (j = ptr1[i], k = ptr2[i]; j < ptr1[i + 1] && k < ptr2[i + 1];) {
                if (ind1[j] == ind2[k]) {
                    if (compVals && da_abs (val1[j] - val2[k]) > eps) {
                        COMPERRPRINT (i + 1, ind1[j] + 1, val1[j], i + 1,
                                ind2[k] + 1, val2[k], rc, fr);
                        ndiff++;
                    }
                    j++;
                    k++;
                } else if (ind1[j] > ind2[k]) {
                    COMPERRPRINT_A (i + 1, ind2[k] + 1, val2[k], rc, fr);
                    k++;
                    ndiff++;
                } else {
                    COMPERRPRINT_B (i + 1, ind1[j] + 1, val1[j], rc, fr);
                    j++;
                    ndiff++;
                }
            }
            for (; j < ptr1[i + 1]; j++) {
                COMPERRPRINT_B (i + 1, ind1[j] + 1, val1[j], rc, fr);
                ndiff++;
            }
            for (; k < ptr2[i + 1]; k++) {
                COMPERRPRINT_A (i + 1, ind2[k] + 1, val2[k], rc, fr);
                ndiff++;
            }
            if (fr)
                printf ("\n");
        }
        for (l = i; i < a->nrows; i++) {
            for (fr = 0, j = ptr1[i]; j < ptr1[i + 1]; j++) {
                COMPERRPRINT_B (i + 1, ind1[j] + 1, val1[j], rc, fr);
                ndiff++;
            }
            if (fr)
                printf ("\n");
        }
        i = l;
        for (; i < b->nrows; i++) {
            for (fr = 0, k = ptr2[i]; k < ptr2[i + 1]; k++) {
                COMPERRPRINT_A (i + 1, ind2[k] + 1, val2[k], rc, fr);
                ndiff++;
            }
            if (fr)
                printf ("\n");
        }

    }

    printf ("Overall, %zu differences were encountered between A and B.\n",
            ndiff);

    da_csr_FreeAll (&a, &b, LTERM);
}


/**
 * Compute 1/sqrt(x) on avg 4x faster than libm
 * http://eggroll.unbsj.ca/rsqrt/rsqrt.pdf
 */
float rsqrt32 (float number)
{
    uint32_t i;
    float x2, y;
    x2 = number * 0.5F;
    y = number;
    i = *(uint32_t *) & y;
    i = 0x5f375a86 - (i >> 1);
    y = *(float *) &i;
    y = y * (1.5F - (x2 * y * y));
    return y;
}

double rsqrt64 (double number)
{
    uint64_t i;
    double x2, y;
    x2 = number * 0.5;
    y = number;
    i = *(uint64_t *) & y;
    i = 0x5fe6eb50c7b537a9 - (i >> 1);
    y = *(double *) &i;
    y = y * (1.5 - (x2 * y * y));
    return y;
}


/**
 * Print compile choices for L2AP
 */
void printCompileChoices ()
{
    printf ("types of pruning: ");
    int i = 0;
#ifdef L2PS
    if (i++)
        printf ("& ");
    printf ("l2ps ");
#elif defined(LENPS)
    if (i++)
        printf ("& ");
    printf ("lenps ");
#endif
#if defined(L2PS) && defined(LENPS)
    da_errexit ("Only one of L2PS and LENPS can be used for index reduction.");
#endif

#ifdef SZ1
    if (i++)
        printf ("& ");
    printf ("sz1 ");
#elif defined(SZ3)
    if (i++)
        printf ("& ");
    printf ("sz3 ");
#endif

#ifdef RS3
    if (i++)
        printf ("& ");
    printf ("rs3 ");
#elif defined(RS1)
    if (i++)
        printf ("& ");
    printf ("rs1 ");
#endif
#ifdef RS2
    if (i++)
        printf ("& ");
    printf ("rs2 ");
#elif defined(RS4)
    if (i++)
        printf ("& ");
    printf ("rs4 ");
#endif
#if defined(RS2) && defined(RS4)
    da_errexit
    ("Only one of RS2 and RS4 can be used for candidate generation pruning.");
#endif

#ifdef L2CG
    if (i++)
        printf ("& ");
    printf ("l2cg ");
#elif defined(LENCG)
    if (i++)
        printf ("& ");
    printf ("lencg ");
#endif
#if defined(L2CG) && defined(LENCG)
    da_errexit
    ("Only one of L2CG and LENCG can be used for candidate generation pruning.");
#endif

#ifdef PSCV
    if (i++)
        printf ("& ");
    printf ("pscv ");
#endif

#ifdef DP1
    if (i++)
        printf ("& ");
    printf ("dp1 ");
#endif
#ifdef DP2
    if (i++)
        printf ("& ");
    printf ("dp2 ");
#endif
#ifdef DP3
    if (i++)
        printf ("& ");
    printf ("dp3 ");
#endif
#ifdef DP4
    if (i++)
        printf ("& ");
    printf ("dp4 ");
#endif
#ifdef DP5
    if (i++)
        printf ("& ");
    printf ("dp5 ");
#endif
#ifdef DP6
    if (i++)
        printf ("& ");
    printf ("dp6 ");
#endif
#ifdef DP7
    if (i++)
        printf ("& ");
    printf ("dp7 ");
#endif
#ifdef DP8
    if (i++)
        printf ("& ");
    printf ("dp8 ");
#endif

#ifdef L2CV
    if (i++)
        printf ("& ");
    printf ("l2cv ");
#elif defined(LENCV)
    if (i++)
        printf ("& ");
    printf ("lencv ");
#endif
#if defined(L2CV) && defined(LENCV)
    da_errexit
    ("Only one of L2CV and LENCV can be used for candidate verification pruning.");
#endif


    printf ("\n");

}


void da_free_sims (da_sims_t ** s)
{
    if ((*s)->sims)
        da_free ((void **) &(*s)->sims, LTERM);
    if ((*s)->ptr)
        da_free ((void **) &(*s)->ptr, LTERM);
    da_free ((void **) s, LTERM);
}
