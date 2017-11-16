/*!
 \file  struct.h
 \brief This file contains structures needed in the program

 \author David C. Anastasiu
 */
#ifndef _L2AP_STRUCT_H_
#define _L2AP_STRUCT_H_

#include "includes.h"

/********************************************************************/
/*! Generator for da_??KeyVal_t data structure */
/********************************************************************/
#define DA_MKKEYVALUE_T(NAME, KEYTYPE, VALTYPE) \
typedef struct {\
  KEYTYPE key;\
  VALTYPE val;\
} NAME;\

/* The actual KeyVal data structures */
	// ptr_t key, various values
	DA_MKKEYVALUE_T(da_ppkv_t,   ptr_t, ptr_t)
	DA_MKKEYVALUE_T(da_pikv_t,   ptr_t, idx_t)
	DA_MKKEYVALUE_T(da_pvkv_t,   ptr_t, val_t)
	DA_MKKEYVALUE_T(da_pckv_t,   ptr_t, char)
	DA_MKKEYVALUE_T(da_pi32kv_t, ptr_t, int32_t)
	DA_MKKEYVALUE_T(da_pi64kv_t, ptr_t, int64_t)
	DA_MKKEYVALUE_T(da_pzkv_t,   ptr_t, ssize_t)
	DA_MKKEYVALUE_T(da_pfkv_t,   ptr_t, float)
	DA_MKKEYVALUE_T(da_pdkv_t,   ptr_t, double)
	DA_MKKEYVALUE_T(da_pskv_t,   ptr_t, char *)
	// idx_t key, various values
	DA_MKKEYVALUE_T(da_ipkv_t,   idx_t, ptr_t)
	DA_MKKEYVALUE_T(da_iikv_t,   idx_t, idx_t)
	DA_MKKEYVALUE_T(da_ivkv_t,   idx_t, val_t)
	DA_MKKEYVALUE_T(da_ickv_t,   idx_t, char)
	DA_MKKEYVALUE_T(da_ii32kv_t, idx_t, int32_t)
	DA_MKKEYVALUE_T(da_ii64kv_t, idx_t, int64_t)
	DA_MKKEYVALUE_T(da_izkv_t,   idx_t, ssize_t)
	DA_MKKEYVALUE_T(da_ifkv_t,   idx_t, float)
	DA_MKKEYVALUE_T(da_idkv_t,   idx_t, double)
	DA_MKKEYVALUE_T(da_iskv_t,   idx_t, char *)
	// size_t key, various values
	DA_MKKEYVALUE_T(da_upkv_t,   size_t, ptr_t)
	DA_MKKEYVALUE_T(da_uikv_t,   size_t, idx_t)
	DA_MKKEYVALUE_T(da_uvkv_t,   size_t, val_t)
	DA_MKKEYVALUE_T(da_uckv_t,   size_t, char)
	DA_MKKEYVALUE_T(da_ui32kv_t, size_t, int32_t)
	DA_MKKEYVALUE_T(da_ui64kv_t, size_t, int64_t)
	DA_MKKEYVALUE_T(da_uzkv_t,   size_t, ssize_t)
	DA_MKKEYVALUE_T(da_ufkv_t,   size_t, float)
	DA_MKKEYVALUE_T(da_udkv_t,   size_t, double)
	DA_MKKEYVALUE_T(da_uskv_t,   size_t, char *)

	// float-float, for use in uniform distributions
	DA_MKKEYVALUE_T(da_ffkv_t,   float, float)


/**
 * A similarity value
 */
typedef struct da_sim_t {
	idx_t i;
	idx_t j;
	float val;
} da_sim_t;

/**
 * A list of similarity values
 */
typedef struct da_sims_t {
	idx_t nsims;
	idx_t size;
	da_sim_t *sims;
} da_sims_t;


/*-------------------------------------------------------------
 * The following data structure stores a sparse CSR format
 *-------------------------------------------------------------*/
typedef struct da_csr_t {
	idx_t nrows, ncols;
	ptr_t *rowptr, *colptr;
	idx_t *rowind, *colind;
	idx_t *rowids, *colids;
	val_t *rowval, *colval;
	val_t *rnorms, *cnorms;
	val_t *rsums, *csums;
	val_t *rsizes, *csizes;
	val_t *rvols, *cvols;
	val_t *rwgts, *cwgts;
} da_csr_t;


/*************************************************************************/
/*! This data structure stores the various variables that make up the 
 * overall state of the system. */
/*************************************************************************/
typedef struct {
	int32_t verbosity;            /* The reporting verbosity level */
    char mode;                    /* What algorithm to execute */
    char sim;                     /* What similarity function to use */
    val_t simT;                   /* Similarity threshold */
    val_t eps;                    /* Adjusted similarity threshold for Tanimoto */
	val_t simT2;                  /* Similarity threshold squared */
    val_t b;                      /* b= .5 + .5/eps */
    val_t b2;                     /* b2= b^2 */
    val_t h;                      /* h= |u|/|v| + |v|/|u| */
    val_t a;                      /* a= (b + sqrt(b^2 - 1)) */
	char scale;                   /* Whether to scale data (by IDF) during pre-processing. */
	char norm;                    /* Whether to normalize the data during pre-processing & what norm to use (l1 vs. l2). */
	char prunerows;               /* Whether to prune rows too long or too short */
	int32_t prminlen;             /* Minimum row length */
	int32_t prmaxlen;             /* Maximum row length */
	char prunecols;               /* Whether to prune cols too long or too short */
	int32_t pcminlen;             /* Minimum col length */
	int32_t pcmaxlen;             /* Maximum col length */
	char compactcols;             /* Compact columns in the matrix */
	char compactrows;             /* Compact rows in the matrix */
	float fldelta;                /* Float delta, for testing matrix value equality. */
	int32_t seed;                 /* Seed for the random number generator */
	char nim;                     /* Whether neighbors should be stored in memory */
	idx_t nqrows;
	idx_t ndrows;
	idx_t ninnzs;

	char fmtRead;                 /* What format the data is stored as: e.g. DA_FMT_BINROW, DA_FMT_CLUTO, DA_FMT_CSR.*/
	char readVals;
	char readNum;
	char fmtWrite;                /* What format the data should be written in: e.g. DA_FMT_BINROW or DA_FMT_CLUTO.*/
	char writeVals;
	char writeNum;

	char *iFile;                  /* The filestem of the input data CSR matrix file. */
	char *oFile;                  /* The filestem of the output file. */
	char *dataset;                /* Dataset name, extracted from filename */
	char *filename;               /* temp space for creating output file names */
	FILE *fpout;                  /* file pointer used for output */
	da_sims_t *neighcache;        /* Neighbor output cache */
	da_csr_t  *docs;              /* Documents structure */
	da_csr_t  *neighbors;         /* Neighbors structure */

	/* internal vars */
	ptr_t nghnnz;                 // current max nnz in the neighbors array
	ptr_t nghinccnt;              // number of times we've increased neighbors
	idx_t progressInd;            // progress indicator chunk
	ssize_t indexSize;            // number of values in the dynamically built inverse index
	ssize_t nSimPairs;            // number of similarity pairs
	ssize_t nDotProducts;         // number of full dot products computed
	ssize_t nCandidates;          // number of candidates considered
	ssize_t nPruneMinsize;        // number of index values pruned by the minsize bound
	ssize_t nPruneDotP;           // number of candidates pruned by the cheap dotp bound
	ssize_t nPruneDotP2;          // number of candidates pruned by the positional dotp bound
	ssize_t nPrunePscore;         // number of candidates pruned by pscore bound
	ssize_t nPruneLength;         // number of candidates pruned by length bound during candidate generation
	ssize_t nPruneLength2;        // number of candidates pruned by length bound during candidate verification
    ssize_t nPruneTanLen;         // number of objects ignored due to original vector length
    ssize_t nPruneTanLenCV;       // number of candidates pruned when checking Tanimoto candidate length in CV
    ssize_t nCountRS;             // number of objects where the Cosine RS bound was effective
    ssize_t nCountTanRS;          // number of objects where the Tanomoto RS bound was effective

	/* bayeslsh related */
	float minMathesT;  			  // threshold for matching, r = c2r(simT)
	int nHashBits;			      // nSketches * NBITS_SKETCH;
	uint32_t *sketches;		      // sketches
	int32_t *minMatches;		  // minimum number of hashes that should be observed to meet simT
	ssize_t *numPruned;			  // number pruned by each sketch limit

	/* row and col permutations */
	idx_t *rperm;                 // row permutation
	idx_t *cperm;                 // col permutation

	/* timers */
	double timer_global;
	double timer_1;
	double timer_2;
	double timer_3;
	double timer_4;
	double timer_5;
	double timer_6;
	double timer_7;
	double timer_8;
	double timer_9;
} params_t;




/**
 * Useful memory for a search task
 */
typedef struct da_tsearch_t {
    idx_t nrows, ncols, nqrows, ndrows;
    ptr_t ninnzs;
    val_t eps;  /* original similarity threshold */
    val_t simT; /* potentially adjusted similarity threshold (for Tan) */
    val_t isimT; /* inverse of simT */
    val_t alpha;
    accum_t *accum;
    val_t *hashval;
    val_t *hashsum;
    val_t *hashlen;
    val_t *hashmax;
    idx_t *hashsz;
    ptr_t *start;
    idx_t *cands;
    da_csr_t *nbrs;
    size_t nbrsz;  // max number of rows / alternatively, max list nnzs
    size_t nbrnnz; // max number of nnzs / alternatively, current list size

    ssize_t nSimPairs;            // number of similarity pairs
    ssize_t nDotProducts;         // number of full dot products computed
    ssize_t nDotProducts1;        // number of full dot products computed in initial construction
    ssize_t nDotProducts2;        // number of full dot products computed in enhancement
    ssize_t nCandidates;          // number of candidates considered
    ssize_t nPruneMinsize;        // number of index values pruned by the minsize bound
    ssize_t nPruneDotP;           // number of candidates pruned by the cheap dotp bound
    ssize_t nPruneDotP2;          // number of candidates pruned by the positional dotp bound
    ssize_t nPrunePscore;         // number of candidates pruned by pscore bound
    ssize_t nPruneLength;         // number of candidates pruned by length bound during candidate generation
    ssize_t nPruneLength2;        // number of candidates pruned by length bound during candidate verification
    ssize_t nPruneTanLen;         // number of non-zeros pruned when initially checking Tanimoto length (0-sim candidates included)
    ssize_t nPruneTanLenCV;       // number of candidates pruned when checking Tanimoto candidate length in CV

    double timer_1;
    double timer_2;
    double timer_3;
    double timer_4;
    double timer_5;
    double timer_6;
    double timer_7;
    double timer_8;
    double timer_9;
} da_tsearch_t;


/***********************************************************************************/
/*! This data structure stores an index for the pl2ap problem. */
/***********************************************************************************/
typedef struct {
    size_t sz;                    // size of index - ptr length is sz+1
    val_t minlen;                 // minimum length of a vector in this index
    val_t maxlen;                 // maximum length of a vector in this index
    ptr_t *ptr;                   // indices in ind/val/len for start of each list
    idx_t *ind;                   // array of ids (list docs)
    val_t *val;                   // array of vals (list values associated with index doc)
    val_t *len;                   // array of lengths associated with vals
    val_t *sum;                   // array of prefix sums with vals
} da_idx_t;


/*************************************************************************/
/*! This data structure stores an inverted index for the all-pairs problem. */
/*************************************************************************/
typedef struct {
	da_ivkv_t *data;              // array of kvf structures (list values)
	da_ivkv_t **starts;			  // pointers in data for start of each list
	da_ivkv_t **ends;             // pointers in data for end of each list
	val_t *accum;                 // accumulator for a given doc
} da_invIdxKv_t;

/*************************************************************************/
/*! This data structure stores an inverted index for the all-pairs problem. */
/*************************************************************************/
typedef struct {
	idx_t nrows;                  // size of accum and cands arrays
	idx_t ncols;                  // size of pointer arrays (starts & ends)
	idx_t nnz;                    // size of index (max num of elements in ids and vals)
	idx_t *ids;                   // array of ids (list docs)
	val_t *vals;                  // array of vals (list values associated with index doc)
	ptr_t *starts;			      // indices in ids/val for start of each list
	ptr_t *ends;                  // indices in ids/val for end of each list
	idx_t *cands;                 // candidates array
	accum_t *accum;               // accumulator for a given doc
} da_invIdx_t;


/*************************************************************************/
/*! This data structure stores an inverted index for the mmjoin all-pairs problem. */
/*************************************************************************/
typedef struct {
	idx_t nrows;                  // size of accum and cands arrays
	idx_t ncols;                  // size of pointer arrays (starts & ends)
	idx_t nnz;                    // size of index (max num of elements in ids and vals)
	idx_t *ids;                   // array of ids (list docs)
	val_t *vals;                  // array of vals (list values associated with index doc)
	ptr_t *starts;			      // indices in ids/val for start of each list
	ptr_t *ends;                  // indices in ids/val for end of each list
	idx_t *cands;                 // candidates array
	val_t *lens;                  // array of sufix lengths associated with vals
	accum_t *accum;               // accumulator for a given doc
	char *mark;
} da_invIdxJ_t;


#endif 
