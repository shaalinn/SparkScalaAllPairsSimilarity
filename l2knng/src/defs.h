/*!
 \file  defs.h
 \brief This file contains various constant and parameter definitions

 \author David C. Anastasiu
 */
#ifndef _L2KNNG_DEFS_H_
#define _L2KNNG_DEFS_H_

/* Versions */
#define PROGRAM_NAME        "knng"
#define VER_MAJOR           0
#define VER_MINOR           1
#define VER_SUBMINOR        0
#define VER_COMMENT         "initial public release"

/** General parameter definitions **/
//#define EXTRACOUNTS  // display extra counts at the end of program execution
//#define EXTRATIMES   // display times for sub-sections of the program execution

#define PRUNENGBS         /* prune neighbor's neighbors below threshold */

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
#define L2PS   // index reduction based on l-2 norm
#define RS4    // use the l2-norm bound to reduce candidate pool during C.G.
#define L2CG   // check accum + ||x_p|| * ||y_p|| during candidate generation
#define L2CV   // check accum + ||x_p|| * ||y_p|| during candidate verification
#define PSCV   // store pscore - last term and use during candidate pruning
#define DP5    // Bayardo's AP bound: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x|,|\vec y'|) \times \rmax_x \times \rmax_y'
#define DP6    // DP5 + prefix max weights: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x|,|\vec y_p'|) \times \rmax_{x_p}' \times \rmax_{y_p}'

/* initial number of neighbors per row to allocate in neighborhood array */
#ifndef NINITNEIGHBORS
	#define NINITNEIGHBORS 10
#endif
/* size of neighbors cache: NSIMSBUFFER * (2 idx + 1 val) */
#ifndef NSIMSBUFFER
	#define NSIMSBUFFER 16384
#endif
#define LENBNDEPS            1e-4         /* epsilon for length bound 2 */

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

//#define IDXCOMPARE       /* compare against IDXJOIN for the splits created by l2knn */
#define IDXLASTSPLIT 0.1  /* execute kidxjoin on the last split or when less than x % of the rows remain */

/* Command-line option codes */
#define CMD_K                   22
#define CMD_MU                  23
#define CMD_STEP                24
#define CMD_NIM                 25
#define CMD_NBL                 27
#define CMD_ALPHA               28
#define CMD_NENHANCE            29
#define CMD_FMT_WRITE           31
#define CMD_FMT_WRITE_NUM       32
#define CMD_WRITE_VALS          33
#define CMD_WRITE_NUMBERING     34
#define CMD_FMT_READ            35
#define CMD_FMT_READ_NUM        36
#define CMD_READ_VALS           37
#define CMD_READ_NUMBERING      38
#define CMD_VERIFY              39
#define CMD_SEED                40
#define CMD_SCALE               51
#define CMD_NORM                52
#define CMD_COMPACT_COLS        53
#define CMD_COMPACT_ROWS        54
#define CMD_SORT_COLS           55
#define CMD_PR                  56
#define CMD_PRMINLEN            57
#define CMD_PRMAXLEN            58
#define CMD_PC                  59
#define CMD_PCMINLEN            60
#define CMD_PCMAXLEN            61
#define CMD_NQROWS              62
#define CMD_NDROWS              63
#define CMD_IGNORE_SIZE         65
#define CMD_PERM_COLS           68
#define CMD_STATS               70
#define CMD_CR                  80
#define CMD_FLDELTA             90
#define CMD_CMPIND              91
#define CMD_CMPVAL              92
#define CMD_NTHREADS            100
#define CMD_VERBOSITY           105
#define CMD_VERSION             109
#define CMD_HELP                110

/* whether to use internal random functions */
#ifndef USE_DARAND
    #define USE_DARAND         1
#endif

/* signal end of list of pointers */
#ifndef LTERM
    #define LTERM  (void **) 0
#endif

/* mopt_t types */
#define DA_MOPT_MARK            1
#define DA_MOPT_CORE            2
#define DA_MOPT_HEAP            3

/* Execution modes */
#define MODE_TESTEQUAL          99  /* Test whether two matrices contain the same values */
#define MODE_IO                 98  /* Transform a matrix from some format into another */
#define MODE_INFO               97  /* Find information about a matrix */
#define MODE_RECALL             96  /* Compute recall given true solution */
#define MODE_IDXJOIN            1   /* IdxJoin */
#define MODE_L2AP               4   /* L2-norm AllPairs */
#define MODE_KIDXJOIN           5   /* IdxJoin solution for K-Nearest Neighbor Graph */
#define MODE_KL2AP              6   /* L2AP solution for K-Nearest Neighbor Graph */
#define MODE_L2KNN              7   /* Efficient K-Nearest Neighbor Graph using L2-Norm bounds */
#define MODE_L2KNN_A            8   /* Approximate solution only for L2KNN */
#define MODE_PKIDXJOIN          9   /* Parallel IdxJoin solution for K-Nearest Neighbor Graph */
#define MODE_MSC               10   /* Maxscore IR solution */
#define MODE_BMM               11   /* Block-Max Maxscore with variable block size */
#define MODE_BMMC              12   /* Block-Max Maxscore with variable block size and compression */
#define MODE_GF                13   /* Greedy Filtering */

/* Permutation modes */
#define PERM_NONE               1   /* do not permute columns before computing */
#define PERM_L2AP               2   /* permute columns as L2AP would do it - decreasing col nnz order */
#define PERM_L2KNN              3   /* permute columns as L2KNN would do it - increasing col nnz order */


/* Similarity type */
#define DA_SIM_COS        1     /* Cosine similarity */
#define DA_SIM_JAC        2     /* Jaccard index/coefficient */
#define DA_SIM_MIN        3     /* Minimum distance/similarity */
#define DA_SIM_AMIN       4     /* Assymetric MIN similarity */
#define DA_SIM_EUC        55    /* Euclidean distance/similarity */
#define DA_SIM_DICE       56    /* Sorensen-Dice coefficient */
#define DA_SIM_OVER       57    /* Overlap distance */
#define DA_SIM_MAN        58    /* Manhattan distance */
#define DA_SIM_TAN        59    /* Tanimoto coefficient */

/* CSR structure components */
#define DA_ROW                  1   /* row-based structure */
#define DA_COL                  2   /* col-based structure */
#define DA_ROWCOL               3   /* both row and col-based */

/* sorting types */
#define DA_SORT_I         1    /* sort in increasing order */
#define DA_SORT_D         2    /* sort in decreasing order */

/* scaling types */
#define DA_SCALE_MAXTF    1    /* TF' = .5 + .5*TF/MAX(TF) */
#define DA_SCALE_MAXTF2   10   /* TF' = .1 + .9*TF/MAX(TF) */
#define DA_SCALE_SQRT     2    /* TF' = .1+SQRT(TF) */
#define DA_SCALE_POW25    3    /* TF' = .1+POW(TF,.25) */
#define DA_SCALE_POW65    4    /* TF' = .1+POW(TF,.65) */
#define DA_SCALE_POW75    5    /* TF' = .1+POW(TF,.75) */
#define DA_SCALE_POW85    6    /* TF' = .1+POW(TF,.85) */
#define DA_SCALE_LOG      7    /* TF' = 1+log_2(TF) */
#define DA_SCALE_IDF      8    /* TF' = TF*IDF */
#define DA_SCALE_IDF2     9    /* TF' = TF*IDF ?? */

/* CSR input formats */
#define DA_FMT_CSR          2
#define DA_FMT_METIS        3
#define DA_FMT_CLUTO        1
#define DA_FMT_BINROW       4
#define DA_FMT_BINCOL       5
#define DA_FMT_IJV          6
#define DA_FMT_BIJV         7
#define DA_FMT_SMAT         50  /* matlab sparse matrix */
#define DA_FMT_BINAP        51  /* satubin - binary format used by allPairs */
#define DA_FMT_BINAPB       52  /* satubin - binary format used by allPairs and PPJoin programs for binary (non-real) input */
#endif
  
 
