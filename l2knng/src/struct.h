/*!
 \file  struct.h
 \brief This file contains structures needed in the program

 \author David C. Anastasiu
 */
#ifndef _L2KNNG_STRUCT_H_
#define _L2KNNG_STRUCT_H_

#include "includes.h"

/********************************************************************/
/*! Generator for da_??KeyVal_t data structure                      */
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
DA_MKKEYVALUE_T(da_iakv_t,   idx_t, accum_t)
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
DA_MKKEYVALUE_T(da_uakv_t,   size_t, accum_t)
DA_MKKEYVALUE_T(da_uckv_t,   size_t, char)
DA_MKKEYVALUE_T(da_ui32kv_t, size_t, int32_t)
DA_MKKEYVALUE_T(da_ui64kv_t, size_t, int64_t)
DA_MKKEYVALUE_T(da_uzkv_t,   size_t, ssize_t)
DA_MKKEYVALUE_T(da_ufkv_t,   size_t, float)
DA_MKKEYVALUE_T(da_udkv_t,   size_t, double)
DA_MKKEYVALUE_T(da_uskv_t,   size_t, char *)


DA_MKKEYVALUE_T(da_uukv_t,   uint, uint)

// float-float, for use in uniform distributions
DA_MKKEYVALUE_T(da_ffkv_t,   float, float)


#define DA_MK2KEYVALUE_T(NAME, KEYTYPE1, KEYTYPE2, VALTYPE) \
typedef struct {\
  KEYTYPE1 key1;\
  KEYTYPE2 key2;\
  VALTYPE val;\
} NAME;\

// iiakv
DA_MK2KEYVALUE_T(da_iiakv_t,  idx_t, idx_t, accum_t)

// Options
/*-------------------------------------------------------------
 * The following data structure implements a string-2-int mapping
 * table used for parsing command-line options
 *-------------------------------------------------------------*/
typedef struct da_StringMap_t {
    char *name;
    int id;
} da_StringMap_t;

extern const da_StringMap_t mode_options[];
extern const da_StringMap_t sim_options[];
extern const da_StringMap_t perm_options[];
extern const da_StringMap_t fmt_options[];
extern const da_StringMap_t scale_options[];


// Some structures from GKlib

/*-------------------------------------------------------------
 * The following data structure stores stores a string as a
 * pair of its allocated buffer and the buffer itself.
 *-------------------------------------------------------------*/
typedef struct da_str_t {
    size_t len;
    char *buf;
} da_str_t;



/*************************************************************************/
/*! The following data structure stores information about a memory
        allocation operation that can either be served from DAmcore_t or by
        a DAmalloc if not sufficient workspace memory is available.      */
/*************************************************************************/
typedef struct da_mop_t {
    int type;
    ssize_t nbytes;
    void *ptr;
} da_mop_t;


/**************************************************************************/
/*! The following structure stores information used for memory allocation */
/**************************************************************************/
typedef struct da_mcore_t {
    /* Workspace information */
    size_t coresize;     /*!< The amount of core memory that has been allocated */
    size_t corecpos;     /*!< Index of the first free location in core */
    void *core;          /*!< Pointer to the core itself */

    /* These are for implementing a stack-based allocation scheme using both
         core and also dynamically allocated memory */
    size_t nmops;         /*!< The number of maop_t entries that have been allocated */
    size_t cmop;          /*!< Index of the first free location in maops */
    da_mop_t *mops;       /*!< The array recording the maop_t operations */

    /* These are for keeping various statistics for wspacemalloc */
    size_t num_callocs;   /*!< The number of core mallocs */
    size_t num_hallocs;   /*!< The number of heap mallocs */
    size_t size_callocs;  /*!< The total # of bytes in core mallocs */
    size_t size_hallocs;  /*!< The total # of bytes in heap mallocs */
    size_t cur_callocs;   /*!< The current # of bytes in core mallocs */
    size_t cur_hallocs;   /*!< The current # of bytes in heap mallocs */
    size_t max_callocs;   /*!< The maximum # of bytes in core mallocs at any given time */
    size_t max_hallocs;   /*!< The maximum # of bytes in heap mallocs at any given time */

} da_mcore_t;



/********************************************************************/
/*! Generator for da_?pq_t data structure                           */
/********************************************************************/
#define DA_MKPQUEUE_T(NAME, KVTYPE)\
typedef struct {\
  size_t nnodes; /* number of items currently in the queue */ \
  size_t maxsize;/* the max size of the queue (additional items dropped) */\
  size_t maxnodes; /* max index of an item stored in the heap */\
\
  /* Heap version of the data structure */ \
  KVTYPE   *heap;  /* array of size maxsize */ \
  ptr_t *locator;  /* array of size maxnodes */ \
} NAME;\

DA_MKPQUEUE_T(da_ivpq_t,    da_ivkv_t)
DA_MKPQUEUE_T(da_iapq_t,    da_iakv_t)

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
	ptr_t *ptr;
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
 * overall state of the system.                                          */
/*************************************************************************/
typedef struct {
	int32_t verbosity;            /* The reporting verbosity level */
	char mode;                    /* What algorithm to execute */
	float simT;                   /* Similarity threshold */
	double simT2;                 /* Similarity threshold squared */
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
	char sortcols;                /* Sort column indices */
	char cmpind;                  /* compare indices in mode testeq */
	char cmpval;                  /* compare values in mode testeq */
	char permcol;                 /* how to permute columns for kij and pkij modes before computing */
	char stats;                   /* Display additional statistics for the matrix in info mode. */
	float fldelta;                /* Float delta, for testing matrix value equality. */
	int32_t seed;                 /* Seed for the random number generator */
	char nim;                     /* Whether neighbors should be stored in memory */
	char cr;                      /* Use checkpoint-restart on enabled methods */
	int32_t nqrows;                /* Size of query block */
	int32_t ndrows;                /* Size of db block - block searching against */
	int32_t sim;                   /* Which similarity to use in computations, e.g., cos, jac. */
	int32_t k;                     /* k in K-NN */
    float alpha;                  /* size of initial approximate neighborhood as a multiple of k */
    float mu;                     /* size of initial approximate neighborhood. mu = k*alpha */
    int32_t pmerge;                /* number of inv index columns to merge at once */
    int32_t nenhance;              /* number of times to try enhance the approximate neighborhood */
    int32_t nthreads;              /* number of threads for OMP parallel execution */
    int32_t nbl;                   /* number of blocks to split into when processing l2knn */
    float step;                   /* step size for knngl2ap */

	char fmtRead;                 /* What format the data is stored as: e.g. DA_FMT_BINROW, DA_FMT_CLUTO, DA_FMT_CSR.*/
	char readVals;
	char readNum;
	char fmtWrite;                /* What format the data should be written in: e.g. DA_FMT_BINROW or DA_FMT_CLUTO.*/
	char writeVals;
	char writeNum;

	char *iFile;                  /* The filestem of the input data CSR matrix file. */
    char *oFile;                  /* The filestem of the output file. */
    char *vFile;                  /* The filestem of the verification file. */
	char *dataset;                /* Dataset name, extracted from filename */
	char *filename;               /* temp space for creating output file names */
    char *crfname;                /* filename for checkpoint restart */
    char *crfname2;               /* filename for checkpoint restart */
	FILE *fpout;                  /* file pointer used for output */
	da_sims_t *neighcache;        /* Neighbor output cache */
    da_csr_t  *docs;              /* Documents structure */
    char *querydocs;              /* Which of the docs we should search neighbors for */
	da_csr_t  *neighbors;         /* Neighbors structure */
	da_iapq_t** knng;             /* K-Nearest Neighbor graph result */

	/* internal vars */
	ptr_t nghnnz;                 // current max nnz in the neighbors array
	ptr_t nghinccnt;              // number of times we've increased neighbors
	idx_t progressInd;            // progress indicator chunk
	ssize_t indexSize;            // number of values in the dynamically built inverse index
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
	val_t   accumPriorPrune;      // accumulated values prior to pruning

	/* row and col permutations */
	idx_t *rperm;                 // row permutation
	idx_t *cperm;                 // col permutation
	idx_t *rsizes;                // row sizes
	idx_t *csizes;                // col sizes
	val_t *rwgts;                 // row max weights
	val_t *cwgts;                 // col max weights

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


/****************************************************************************/
/*! This data structure stores an inverted index for the all-pairs problem. */
/****************************************************************************/
typedef struct {
	da_ivkv_t *data;              // array of kvf structures (list values)
	da_ivkv_t **starts;			  // pointers in data for start of each list
	da_ivkv_t **ends;             // pointers in data for end of each list
	val_t *accum;                 // accumulator for a given doc
} da_invIdxKv_t;

/****************************************************************************/
/*! This data structure stores an inverted index for the all-pairs problem. */
/****************************************************************************/
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


/***********************************************************************************/
/*! This data structure stores an inverted index for the mmjoin all-pairs problem. */
/***********************************************************************************/
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
} da_invIdxJ_t;


/*********************************************************************************************/
/*! This data structure is an extended inverted index - stores input & space for sim search. */
/*********************************************************************************************/
typedef struct {
    da_csr_t *docs;               // input or partial input
    idx_t nrows;                  // size of accum and cands arrays
    idx_t ncols;                  // size of pointer arrays (starts & ends)
    idx_t nnz;                    // size of index (max num of elements in ids and vals)
    idx_t sid;                    // slice id/index id
    char *tag;                   // tag denoting whether row id is no longer being considered in the search
    idx_t *ids;                   // translation of ids to global space
    val_t *lens;                  // array of row suffix lengths
    idx_t *idxind;                // array of ids (list docs)
    val_t *idxval;                // array of vals (list values associated with index doc)
    ptr_t *idxptr;                // indices in ids/val for start of each list
    val_t *idxlens;               // array of sufix lengths associated with vals
    idx_t *cands;                 // candidates array
    val_t *accum;                 // accumulator for a given doc
} da_sfind_t;


/**************************************************************************/
/*! This data structure is a pointer into an inverted list - for Maxscore */
/**************************************************************************/
typedef struct {
    idx_t tid;          // the term id
    idx_t qfr;          // query vector frequency for this col
    float idf;         // idf for this column
    float maxscore;    // the max score of the entire list
    size_t size;        // number of values
    idx_t *ids;         // pointer to ids
    idx_t *freqs;       // pointer to frequencies
} da_il_t;

/**************************************************************************************/
/*! This data structure is a pointer into a (compressed) block inverted list - for BMM. */
/**************************************************************************************/
typedef struct {
    idx_t tid;          // the term id
    idx_t qfr;          // query vector frequency for this col
    float idf;         // idf for this column
    float maxscore;    // the max score of the entire list
    uint did;           // current did being considered in list
    float max;         // max score in block currently considered in list
    uint   bxp;         // block size exponent for this list
    idx_t  nbl;         // number of blocks
    idx_t  i;           // current block iterator
    idx_t  j;           // current location in (de-compressed) ids and frequencies
    char   cfl;         // flag specifying whether docs and/or frequencies have been de-compressed for current block
    uint sz;            // size of de-compressed block

    uint  *dids;        // (optionally compressed) block ids for all blocks in list
    uint  *dfreqs;      // (optionally compressed) block frequencies for all blocks in list
    ptr_t *bptr;        // block pointer in bids and bfreqs or compressed data, delimiting blocks. of size (#blocks+1)
    uint  *bsz;         // number of elements in block (only for compressed version)
    uint  *bdmax;       // block max doc id - of size (#blocks)
    float *bsmax;       // block max scores - of size (#blocks)

    uint  *ids;         // de-compressed ids of current block
    uint  *freqs;       // de-compressed frequencies of current block
} da_bil_t;


#endif 
