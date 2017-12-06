/*!
 \file  defs.h
 \brief This file contains various constant and parameter definitions

 \author David C. Anastasiu
 */
#ifndef _COSGRAPH_DEFS_H_
#define _COSGRAPH_DEFS_H_

/* Versions */
#define PROGRAM_NAME        "pl2ap"
#define VER_MAJOR           0
#define VER_MINOR           0
#define VER_SUBMINOR        11
#define VER_COMMENT         "dev release - added pl2nn"

/** General parameter definitions **/
//#define EXTRACOUNTS   // display extra counts at the end of program execution
#define EXTRATIMES    // display times for sub-sections of the program execution
//#define HASHCOUNTS    // count number of hash collisions
//#define ACCUMCOUNTS   // track the percent of accumulator operations performed before pruning

//#define IDX_SIZE      // report index sizes and exit

#define NITERS          20

/* How to handle forward index dot-product in L2AP */
#define L2APDP_HASH      1  /** Create dense version of query vector */
#define L2APDP_MASK      2  /** Hash masked version of query vector -- linear search for collisions */
#define L2APDP_MIX       3  /** For each query row, choose the best option based on length */
#define L2APDP           L2APDP_MIX

#define PL2APRS_STORE    1  /** Pre-compute and store RS scores */
#define PL2APRS_LIVE     2  /** Compute RS scores on the fly */
#define PL2APRS          PL2APRS_STORE

#define PL2AP_DP6        // pruning based on prefix max weights
#define PL2AP_MINSZ      // pruning based on minimum size
//#define PL2AP_SUP        // update start pointers in index
#ifndef PL2AP_MINSZ
#undef PL2AP_SUP
#endif

/** Length of masking hash table array */
#define HTSIZE          (1<<13)
#define HTMASK          (HTSIZE-1)
#define HTUBND          (HTSIZE-1)>>3

#define PRUNENGBS        /* prune neighbor's neighbors below threshold */

#define COLSORTBY_SIZE  1
#define COLSORTBY_MAX   2
#define COLSORTBY_AVG   3
#define COLSORTBY_RAND  4
#define COLSORTBY       COLSORTBY_SIZE

#define COLSORTORDER_I  1
#define COLSORTORDER_D  2

#define COLSORTORDER    COLSORTORDER_I

#define PMERGE_ALL      1
#define PMERGE_TWOS     2


/** parameter definitions for apss **/
#define L2PS        // index reduction based on l-2 norm
#define RS4         // use the l2-norm bound to reduce candidate pool during C.G.
#define RS1         // use Bayardo's max vector bound to reduce candidate pool during C.G.
#define L2CG        // check accum + ||x_p|| * ||y_p|| during candidate generation
#define L2CV        // check accum + ||x_p|| * ||y_p|| during candidate verification
#define PSCV        // store pscore - last term and use during candidate pruning
#define DP5         // Bayardo's AP bound: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x|,|\vec y'|) \times \rmax_x \times \rmax_y'
#define DP6         // DP5 + prefix max weights: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x|,|\vec y_p'|) \times \rmax_{x_p}' \times \rmax_{y_p}'

/* initial number of neighbors per row to allocate in neighborhood array */
#ifndef NINITNEIGHBORS
#define NINITNEIGHBORS 10
#endif
/* size of neighbors cache: NSIMSBUFFER * (2 idx + 1 val) */
#ifndef NSIMSBUFFER
#define NSIMSBUFFER 16384
#endif
#define LENBNDEPS            1e-4    /* epsilon for length bound 2 */

/* housekeeping for pruning params */
#if defined(L2CG) || defined(L2CV) || defined(RS4)
#define IDXL2
#endif
#if defined(LENCG) || defined(LENCV) || defined(RS2)
#define IDXLEN
#endif
#if defined(SZ1) && !defined(SZ3)
#define SZ3
#endif
#if defined(RS3) && !defined(RS1)
#define RS1
#endif

//#define IDXCOMPARE        /* compare against IDXJOIN for the splits created by l2knn */
#define IDXLASTSPLIT 0.1    /* execute kidxjoin on the last split or when less than x % of the rows remain */

/* Command-line option codes */
#define CMD_EPS           21
#define CMD_K             22
#define CMD_NIM           25
#define CMD_SYMOUT        30
#define CMD_FMT_WRITE     31
#define CMD_FMT_WRITE_NUM 32
#define CMD_WRITE_VALS    33
#define CMD_WRITE_NUM     34
#define CMD_FMT_READ      35
#define CMD_FMT_READ_NUM  36
#define CMD_READ_VALS     37
#define CMD_READ_NUM      38
#define CMD_VERIFY        39
#define CMD_SEED          40
#define CMD_NQROWS        41
#define CMD_NDROWS        42
#define CMD_NINNZ         43
#define CMD_SCALE         51
#define CMD_NORM          52
#define CMD_COMPACT_COLS  53
#define CMD_COMPACT_ROWS  54
#define CMD_SORT_COLS     55
#define CMD_PR            56
#define CMD_PRMINLEN      57
#define CMD_PRMAXLEN      58
#define CMD_PC            59
#define CMD_PCMINLEN      60
#define CMD_PCMAXLEN      61
#define CMD_IGNORE_SIZE   65
#define CMD_STATS         70
#define CMD_CR            80
#define CMD_FLDELTA       90
#define CMD_CMPIND        91
#define CMD_CMPVAL        92
#define CMD_NTHREADS      100
#define CMD_RUN           102
#define CMD_VERBOSITY     105
#define CMD_VERSION       109
#define CMD_HELP          110

/* whether to use internal random functions */
#ifndef USE_DARAND
#define USE_DARAND        1
#endif

/* signal end of list of pointers */
#ifndef LTERM
#define LTERM  (void **) 0
#endif

/* mopt_t types */
#define DA_MOPT_MARK      1
#define DA_MOPT_CORE      2
#define DA_MOPT_HEAP      3

/* Execution modes */
#define MODE_TESTEQUAL    99 /* Test whether two matrices contain the same values */
#define MODE_IO           98 /* Transform a matrix from some format into another */
#define MODE_INFO         97 /* Find information about a matrix */
#define MODE_RECALL       96 /* Compute recall given true solution */
#define MODE_IDXJOIN      1  /* IdxJoin solution for APSS */
#define MODE_PIDXJOIN     2  /* Parallel IdxJoin solution for APSS */
#define MODE_L2AP         3  /* L2-norm AllPairs solution for APSS */
#define MODE_PL2AP        4  /* Parallel L2-norm AllPairs solution for APSS */
#define MODE_L2AP_A       5  /* Approximate L2-norm AllPairs solution for APSS */
#define MODE_PL2AP_A      6  /* Parallel Approximate L2-norm AllPairs solution for APSS */
#define MODE_PIJNN        7  /* Parallel IdxJoin solution for Nearest Neighbor Search */
#define MODE_PL2NN        8  /* Parallel L2-norm Nearest Neighbor search solution given minimum similarity threshold eps */
#define MODE_PL2NN_A      9  /* Parallel Approximate L2-norm Nearest Neighbor search solution given minimum similarity threshold eps */

/* Similarity type */
#define DA_SIM_COS        1     /* Cosine similarity */
#define DA_SIM_JAC        2     /* Jaccard index/coefficient */
#define DA_SIM_MIN        3     /* Minimum distance/similarity */
#define DA_SIM_AMIN       4  /* Assymetric MIN similarity */
#define DA_SIM_EUC        55 /* Euclidean distance/similarity */
#define DA_SIM_DICE       56 /* Sorensen-Dice coefficient */
#define DA_SIM_OVER       57 /* Overlap distance */
#define DA_SIM_MAN        58 /* Manhattan distance */
#define DA_SIM_TAN        59 /* Tanimoto coefficient */

/* CSR structure components */
#define DA_ROW            1     /* row-based structure */
#define DA_COL            2     /* col-based structure */
#define DA_ROWCOL         3     /* both row and col-based */

/* sorting types */
#define DA_SORT_I         1     /* sort in increasing order */
#define DA_SORT_D         2     /* sort in decreasing order */

/* scaling types */
#define DA_SCALE_MAXTF    1     /* TF' = .5 + .5*TF/MAX(TF) */
#define DA_SCALE_MAXTF2   10    /* TF' = .1 + .9*TF/MAX(TF) */
#define DA_SCALE_SQRT     2     /* TF' = .1+SQRT(TF) */
#define DA_SCALE_POW25    3     /* TF' = .1+POW(TF,.25) */
#define DA_SCALE_POW65    4     /* TF' = .1+POW(TF,.65) */
#define DA_SCALE_POW75    5     /* TF' = .1+POW(TF,.75) */
#define DA_SCALE_POW85    6     /* TF' = .1+POW(TF,.85) */
#define DA_SCALE_LOG      7     /* TF' = 1+log_2(TF) */
#define DA_SCALE_IDF      8     /* TF' = TF*IDF */
#define DA_SCALE_IDF2     9     /* TF' = TF*IDF ?? */

/* CSR input formats */
#define DA_FMT_CLUTO      1
#define DA_FMT_CSR        2
#define DA_FMT_METIS      3
#define DA_FMT_BINROW     4
#define DA_FMT_BINCOL     5
#define DA_FMT_IJV        6
#define DA_FMT_BIJV       7
#define DA_FMT_SMAT       50    /* matlab sparse matrix */
#define DA_FMT_BINAP      51    /* satubin - binary format used by allPairs */
#define DA_FMT_BINAPB     52    /* satubin - binary format used by allPairs and PPJoin programs for binary (non-real) input */



#endif
