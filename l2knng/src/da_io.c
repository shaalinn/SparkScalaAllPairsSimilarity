/*!
\file  io.c
\brief Various file I/O functions.

This file contains various functions that perform I/O.

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Aug 2013.

\date Started 4/10/95
\author George
\version\verbatim $Id: io.c 15225 2013-09-25 14:49:08Z karypis $ \endverbatim
 */

#ifdef HAVE_GETLINE
/* Get getline to be defined. */
#define _GNU_SOURCE
#include <stdio.h>
#undef _GNU_SOURCE
#endif

#include "includes.h"

/*************************************************************************
 * This function opens a file
 **************************************************************************/
FILE* da_fopen(const char* const fname, const char* const mode, const char* const msg)
{
	FILE *fp;
	char errmsg[8192];

	fp = fopen(fname, mode);
	if (fp != NULL)
		return fp;

	sprintf(errmsg,"file: %s, mode: %s, [%s]", fname, mode, msg);
	perror(errmsg);
	da_errexit("Failed on da_fopen()\n");

	return NULL;
}


/*************************************************************************
 * This function closes a file
 **************************************************************************/
void da_fclose(FILE* fp)
{
	fclose(fp);
}


/*************************************************************************/
/*! This function is a wrapper around the read() function that ensures 
    that all data is been read, by issueing multiple read requests.
 */
/*************************************************************************/
ssize_t da_read(const int fd, void* vbuf, const size_t count)
{
	char *buf = (char *)vbuf;
	ssize_t rsize, tsize=count;

	do {
		if ((rsize = read(fd, buf, tsize)) == -1)
			return -1;
		buf   += rsize;
		tsize -= rsize;
	} while (tsize > 0);

	return count;
}



/*************************************************************************/
/*! This function is the GKlib implementation of glibc's getline()
    function.
    \returns -1 if the EOF has been reached, otherwise it returns the 
             number of bytes read.
 */
/*************************************************************************/
ssize_t da_getline(char** lineptr, size_t* n, FILE* stream)
{
#ifdef HAVE_GETLINE
	return getline(lineptr, n, stream);
#else
	size_t i;
	int ch;

	if (feof(stream))
		return -1;

	/* Initial memory allocation if *lineptr is NULL */
	if (*lineptr == NULL || *n == 0) {
		*n = 1024;
		*lineptr = da_malloc((*n)*sizeof(char), "da_getline: lineptr");
	}

	/* get into the main loop */
	i = 0;
	while ((ch = getc(stream)) != EOF) {
		(*lineptr)[i++] = (char)ch;

		/* reallocate memory if reached at the end of the buffer. The +1 is for '\0' */
		if (i+1 == *n) {
			*n = 2*(*n);
			*lineptr = da_realloc(*lineptr, (*n)*sizeof(char), "da_getline: lineptr");
		}

		if (ch == '\n')
			break;
	}
	(*lineptr)[i] = '\0';

	return (i == 0 ? -1 : i);
#endif
}


/*************************************************************************/
/*! This function reads the contents of a text file and returns it in the
    form of an array of strings.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
 */
/*************************************************************************/
char** da_readfile(const char* const fname, size_t* r_nlines)
{
	size_t lnlen=0, nlines=0;
	char *line=NULL, **lines=NULL;
	FILE *fpin;

	da_getfilestats(fname, &nlines, NULL, NULL, NULL);
	if (nlines > 0) {
		lines = (char **)da_malloc(nlines*sizeof(char *), "da_readfile: lines");

		fpin = da_fopen(fname, "r", "da_readfile");
		nlines = 0;
		while (da_getline(&line, &lnlen, fpin) != -1) {
			da_strtprune(line, "\n\r");
			lines[nlines++] = da_strdup(line);
		}
		da_fclose(fpin);
	}

	da_free((void **)&line, LTERM);

	if (r_nlines != NULL)
		*r_nlines  = nlines;

	return lines;
}

/*************************************************************************/
/*! This function reads the contents of a file and returns it in the
    form of an array of int32_t.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
 */
/*************************************************************************/
int32_t* da_i32readfile(const char* const fname, size_t* r_nlines)
{
	size_t lnlen=0, nlines=0;
	char *line=NULL;
	int32_t *array=NULL;
	FILE *fpin;

	da_getfilestats(fname, &nlines, NULL, NULL, NULL);
	if (nlines > 0) {
		array = da_i32malloc(nlines, "da_i32readfile: array");

		fpin = da_fopen(fname, "r", "da_readfile");
		nlines = 0;

		while (da_getline(&line, &lnlen, fpin) != -1) {
			sscanf(line, "%"SCNd32, &array[nlines++]);
		}

		da_fclose(fpin);
	}

	da_free((void **)&line, LTERM);

	if (r_nlines != NULL)
		*r_nlines  = nlines;

	return array;
}

/*************************************************************************/
/*! This function reads the contents of a file and returns it in the
    form of an array of int64_t.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
 */
/*************************************************************************/
int64_t* da_i64readfile(const char* const fname, size_t* r_nlines)
{
	size_t lnlen=0, nlines=0;
	char *line=NULL;
	int64_t *array=NULL;
	FILE *fpin;

	da_getfilestats(fname, &nlines, NULL, NULL, NULL);
	if (nlines > 0) {
		array = da_i64malloc(nlines, "da_i64readfile: array");

		fpin = da_fopen(fname, "r", "da_readfile");
		nlines = 0;

		while (da_getline(&line, &lnlen, fpin) != -1) {
			sscanf(line, "%"SCNd64, &array[nlines++]);
		}

		da_fclose(fpin);
	}

	da_free((void **)&line, LTERM);

	if (r_nlines != NULL)
		*r_nlines  = nlines;

	return array;
}

/*************************************************************************/
/*! This function reads the contents of a binary file and returns it in the
    form of an array of int32_t.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
 */
/*************************************************************************/
int32_t* da_i32readfilebin(const char* const fname, size_t* r_nelmnts)
{
	size_t nelmnts;
	ssize_t fsize;
	int32_t *array=NULL;
	FILE *fpin;

	*r_nelmnts = 0;

	fsize = da_getfsize(fname);

	if (fsize == -1) {
		da_errexit( "Failed to fstat(%s).\n", fname);
		return NULL;
	}

	if (fsize%sizeof(int32_t) != 0) {
		da_errexit( "The size [%zd] of the file [%s] is not in multiples of sizeof(int32_t).\n", fsize, fname);
		return NULL;
	}

	nelmnts = fsize/sizeof(int32_t);
	array = da_i32malloc(nelmnts, "da_i32readfilebin: array");

	fpin = da_fopen(fname, "rb", "da_i32readfilebin");

	if (fread(array, sizeof(int32_t), nelmnts, fpin) != nelmnts) {
		da_errexit( "Failed to read the number of words requested. %zd\n", nelmnts);
		da_free((void **)&array, LTERM);
		return NULL;
	}
	da_fclose(fpin);

	*r_nelmnts = nelmnts;

	return array;
}

/*************************************************************************/
/*! This function writes the contents of an array into a binary file.
    \param fname is the name of the file
    \param n the number of elements in the array.
    \param a the array to be written out.
 */
/*************************************************************************/
size_t da_i32writefilebin(const char* const fname, const size_t n, const int32_t* const a)
{
	size_t fsize;
	FILE *fp;

	fp = da_fopen(fname, "wb", "da_writefilebin");

	fsize = fwrite(a, sizeof(int32_t), n, fp);

	da_fclose(fp);

	return fsize;
}

/*************************************************************************/
/*! This function reads the contents of a binary file and returns it in the
    form of an array of int64_t.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
 */
/*************************************************************************/
int64_t* da_i64readfilebin(const char* const fname, size_t* r_nelmnts)
{
	size_t nelmnts;
	ssize_t fsize;
	int64_t *array=NULL;
	FILE *fpin;

	*r_nelmnts = 0;

	fsize = da_getfsize(fname);

	if (fsize == -1) {
		da_errexit( "Failed to fstat(%s).\n", fname);
		return NULL;
	}

	if (fsize%sizeof(int64_t) != 0) {
		da_errexit( "The size of the file is not in multiples of sizeof(int64_t).\n");
		return NULL;
	}

	nelmnts = fsize/sizeof(int64_t);
	array = da_i64malloc(nelmnts, "da_i64readfilebin: array");

	fpin = da_fopen(fname, "rb", "da_i64readfilebin");

	if (fread(array, sizeof(int64_t), nelmnts, fpin) != nelmnts) {
		da_errexit( "Failed to read the number of words requested. %zd\n", nelmnts);
		da_free((void **)&array, LTERM);
		return NULL;
	}
	da_fclose(fpin);

	*r_nelmnts = nelmnts;

	return array;
}

/*************************************************************************/
/*! This function writes the contents of an array into a binary file.
    \param fname is the name of the file
    \param n the number of elements in the array.
    \param a the array to be written out.
 */
/*************************************************************************/
size_t da_i64writefilebin(const char* const fname, const size_t n, const int64_t* const a)
{
	size_t fsize;
	FILE *fp;

	fp = da_fopen(fname, "wb", "da_writefilebin");

	fsize = fwrite(a, sizeof(int64_t), n, fp);

	da_fclose(fp);

	return fsize;
}

/*************************************************************************/
/*! This function reads the contents of a binary file and returns it in the
    form of an array of float.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
 */
/*************************************************************************/
float* da_freadfilebin(const char* const fname, size_t* r_nelmnts)
{
	size_t nelmnts;
	ssize_t fsize;
	float *array=NULL;
	FILE *fpin;

	*r_nelmnts = 0;

	fsize = da_getfsize(fname);

	if (fsize == -1) {
		da_errexit( "Failed to fstat(%s).\n", fname);
		return NULL;
	}

	if (fsize%sizeof(float) != 0) {
		da_errexit( "The size of the file is not in multiples of sizeof(float).\n");
		return NULL;
	}

	nelmnts = fsize/sizeof(float);
	array = da_fmalloc(nelmnts, "da_freadfilebin: array");

	fpin = da_fopen(fname, "rb", "da_freadfilebin");

	if (fread(array, sizeof(float), nelmnts, fpin) != nelmnts) {
		da_errexit( "Failed to read the number of words requested. %zd\n", nelmnts);
		da_free((void **)&array, LTERM);
		return NULL;
	}
	da_fclose(fpin);

	*r_nelmnts = nelmnts;

	return array;
}

/*************************************************************************/
/*! This function writes the contents of an array into a binary file.
    \param fname is the name of the file
    \param n the number of elements in the array.
    \param a the array to be written out.
 */
/*************************************************************************/
size_t da_fwritefilebin(const char* const fname, const size_t n, const float* const a)
{
	size_t fsize;
	FILE *fp;

	fp = da_fopen(fname, "wb", "da_fwritefilebin");

	fsize = fwrite(a, sizeof(float), n, fp);

	da_fclose(fp);

	return fsize;
}

/*************************************************************************/
/*! This function reads the contents of a binary file and returns it in the
    form of an array of double.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
 */
/*************************************************************************/
double* da_dreadfilebin(const char* const fname, size_t* r_nelmnts)
{
	size_t nelmnts;
	ssize_t fsize;
	double *array=NULL;
	FILE *fpin;

	*r_nelmnts = 0;

	fsize = da_getfsize(fname);

	if (fsize == -1) {
		da_errexit( "Failed to fstat(%s).\n", fname);
		return NULL;
	}

	if (fsize%sizeof(double) != 0) {
		da_errexit( "The size of the file is not in multiples of sizeof(double).\n");
		return NULL;
	}

	nelmnts = fsize/sizeof(double);
	array = da_dmalloc(nelmnts, "da_dreadfilebin: array");

	fpin = da_fopen(fname, "rb", "da_dreadfilebin");

	if (fread(array, sizeof(double), nelmnts, fpin) != nelmnts) {
		da_errexit( "Failed to read the number of words requested. %zd\n", nelmnts);
		da_free((void **)&array, LTERM);
		return NULL;
	}
	da_fclose(fpin);

	*r_nelmnts = nelmnts;

	return array;
}

/*************************************************************************/
/*! This function writes the contents of an array into a binary file.
    \param fname is the name of the file
    \param n the number of elements in the array.
    \param a the array to be written out.
 */
/*************************************************************************/
size_t da_dwritefilebin(const char* const fname, const size_t n, const double* const a)
{
	size_t fsize;
	FILE *fp;

	fp = da_fopen(fname, "wb", "da_writefilebin");

	fsize = fwrite(a, sizeof(double), n, fp);

	da_fclose(fp);

	return fsize;
}

