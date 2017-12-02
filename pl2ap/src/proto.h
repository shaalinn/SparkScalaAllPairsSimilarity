/*!
 \file  proto.h
 \brief This file contains function prototypes

 \author David C. Anastasiu
 */
#ifndef _COSGRAPH_PROTO_H_
#define _COSGRAPH_PROTO_H_

#include "includes.h"

#ifdef __cplusplus
extern "C"
{
#endif



/* main.c */
void readInputData (params_t * params);
void preProcessData (params_t * params);
void da_testMatricesEqual (params_t * params);
void da_matrixInfo (params_t * params);
void da_matrixIo (params_t * params);
void simSearchSetup (params_t * params);
void simSearchFinalize (params_t * params);
void freeParams (params_t ** params);



/* l2ap.c */
void l2apFindNeighbors (params_t * params);
void l2apReorderDocs (params_t * params, da_csr_t ** docs, idx_t * rsizes,
		idx_t * csizes, val_t * rwgts, val_t * cwgts,
		idx_t * rperm, idx_t * cperm);
void l2apReorderCols (da_csr_t * docs, idx_t * csizes, idx_t * cperm);
void l2knnReorderCols (da_csr_t * docs, idx_t * csizes, idx_t * cperm);


/* pl2ap.c */
void pl2apFindNeighbors (params_t * params);
float computePercAccum (idx_t qid, idx_t cid, idx_t iic, idx_t fic,
		da_csr_t * docs, ptr_t * endptr);
void pl2apFindMatches (idx_t lrid, idx_t qstart, idx_t dstart, val_t simT,
		da_csr_t * docs, val_t * rs, val_t * lens,
		val_t * prmax, val_t * ps, val_t * rfiwgts,
		val_t * cwgts, val_t rwgt, idx_t * rsizes,
		da_idx_t * iidx, da_idx_t * fidx, da_tsearch_t * ts,
		ptr_t * endptr);
idx_t pl2apProcessCandidates (idx_t rid, idx_t qstart, idx_t dstart,
		val_t simT, size_t sz, val_t rwgt,
		val_t * rfiwgts, val_t * ps, val_t * hashval,
		val_t * hashlen, val_t * hashmax,
		da_idx_t * fidx, idx_t ncands, idx_t * cands,
		accum_t * accum, da_tsearch_t * ts,
		da_csr_t * docs, ptr_t * endptr);
void pl2apFindMatchesMask (idx_t rid, idx_t qstart, idx_t dstart,
		val_t simT, da_csr_t * docs, val_t * rs,
		val_t * lens, val_t * prmax, val_t * ps,
		val_t * rfiwgts, val_t * cwgts, val_t rwgt,
		idx_t * rsizes, da_idx_t * iidx, da_idx_t * fidx,
		da_tsearch_t * ts, ptr_t * endptr);
idx_t pl2apProcessCandidatesMask (idx_t rid, idx_t qstart, idx_t dstart,
		val_t simT, size_t sz, val_t rwgt,
		val_t * rfiwgts, val_t * ps,
		da_idx_t * fidx, size_t hashSize,
		idx_t ncands, idx_t * cands,
		accum_t * accum, da_tsearch_t * ts,
		da_csr_t * docs, ptr_t * endptr);
void pl2apFindIndexSplit (idx_t rid, params_t * params, da_csr_t * docs,
		ptr_t * endptr, val_t * rwgts, val_t * rfiwgts,
		val_t * cwgts, val_t * ps, val_t * rs,
		val_t * lens, val_t * prmax);
void pl2apIndexRows (idx_t start, idx_t nrows, params_t * params,
		da_csr_t * docs, ptr_t * endptr, da_idx_t * iidx,
		da_idx_t * fidx);

/* pl2nn.c */
void pl2nnFindNeighbors (params_t * params);
float computePercAccum (idx_t qid, idx_t cid, idx_t iic, idx_t fic,
		da_csr_t * docs, ptr_t * endptr);
void pl2nnFindMatches (idx_t lrid, idx_t qstart, idx_t dstart, val_t simT,
		da_csr_t * docs, val_t * rs, val_t * lens,
		val_t * prmax, val_t * ps, val_t * rfiwgts,
		val_t * cwgts, val_t rwgt, idx_t * rsizes,
		da_idx_t * iidx, da_idx_t * fidx, da_tsearch_t * ts,
		ptr_t * endptr);
idx_t pl2nnProcessCandidates (idx_t rid, idx_t qstart, idx_t dstart,
		val_t simT, size_t sz, val_t rwgt,
		val_t * rfiwgts, val_t * ps, val_t * hashval,
		val_t * hashlen, val_t * hashmax,
		da_idx_t * fidx, idx_t ncands, idx_t * cands,
		accum_t * accum, da_tsearch_t * ts,
		da_csr_t * docs, ptr_t * endptr);
void pl2nnFindIndexSplit (idx_t rid, params_t * params, da_csr_t * docs,
		ptr_t * endptr, val_t * rwgts, val_t * rfiwgts,
		val_t * cwgts, val_t * ps, val_t * rs,
		val_t * lens, val_t * prmax);
void pl2nnIndexRows (idx_t start, idx_t nrows, params_t * params,
		da_csr_t * docs, ptr_t * endptr, da_idx_t * iidx,
		da_idx_t * fidx);





/* idxjoin.c */
void ijFindNeighbors (params_t * params);
idx_t da_ApssGetSimilarSmallerRows (da_csr_t * mat, idx_t rid,
		char noSelfSim, idx_t nqterms,
		idx_t * qind, val_t * qval,
		char simtype, idx_t nsim, float minsim,
		da_ivkv_t * hits, idx_t * i_marker,
		da_ivkv_t * i_cand, idx_t * ncands);
idx_t da_ApssGetSimilarRows (da_csr_t * mat, idx_t rid, char noSelfSim,
		idx_t nqterms, idx_t * qind, val_t * qval,
		char simtype, idx_t nsim, float minsim,
		da_ivkv_t * hits, idx_t * i_marker,
		da_ivkv_t * i_cand, idx_t * ncands);

/* pidxjoin.c */
void pijFindNeighbors (params_t * params);

/* pijnn.c */
void pijnnFindNeighbors (params_t * params);

/* util.c */
void da_errexit (const char *const f_str, ...);
double da_WClockSeconds (void);
double da_CPUSeconds (void);
char da_getFileFormat (char *file, const char format);
char *da_getDataset (params_t * params);
void da_gWriteVector (char *filename, val_t * vec, idx_t size,
		char *separator);
void da_vWriteVector (char *filename, val_t * vec, idx_t size,
		char *separator);
void da_vWriteMatrix (char *filename, val_t ** mat, idx_t nrows,
		idx_t ncols);
void da_dWriteMatrix (char *filename, double **mat, idx_t nrows,
		idx_t ncols);
void da_iWriteMatrix (char *filename, idx_t ** mat, idx_t nrows,
		idx_t ncols);
void da_iWriteVector (char *filename, idx_t * vec, idx_t size,
		char *separator);
void da_pWriteVector (char *filename, idx_t * vec, idx_t size,
		char *separator);
void da_storeSim (params_t * params, idx_t i, idx_t j, val_t v);
char da_isFmtBinary (char fmt);
void da_printTimer (char *name, double time);
void da_printTimerLong (char *name, double time);
void da_addHigherNeighbors (da_csr_t * mat);
int da_log2 (idx_t a);
int da_ispow2 (idx_t a);
float da_flog2 (float a);
double da_dlog2 (double a);
val_t da_vlog2 (val_t a);
void da_csrCompare (da_csr_t * a, da_csr_t * b, float eps, char compInds,
		char compVals);
float da_fsqrt (float number);
da_csr_t *da_csr_CombineAndPermute (idx_t nrows, idx_t sz,
		da_tsearch_t ** tsearch, idx_t * perm);
da_csr_t *da_csr_CombineAndPermuteLT (idx_t nrows, idx_t sz,
		da_tsearch_t ** tsearch,
		idx_t * perm);
void da_inversePermuteSqMatrixLT (da_csr_t * mat, idx_t * perm, char which);
void da_inversePermuteMatrix (da_csr_t ** matP, idx_t * rowPerm,
		idx_t * colPerm);
void da_inversePermuteRows (da_csr_t * mat, idx_t * rowPerm,
		idx_t * colPerm);
double da_CosineSimilarity (const size_t nind1, const idx_t const *ind1,
		const val_t const *val1, const size_t nind2,
		const idx_t const *ind2,
		const val_t const *val2);
void printCompileChoices ();
void da_free_sims (da_sims_t ** s);

/* cmdline.c */
void cmdline_parse (params_t * ctrl, int argc, char *argv[]);


/***** Library functions ******/

/* omp.c */
#if !defined(_OPENMP)
void omp_set_num_threads (int num_threads);
int omp_get_num_threads (void);
int omp_get_max_threads (void);
int omp_get_thread_num (void);
int omp_get_num_procs (void);
int omp_in_parallel (void);
void omp_set_dynamic (int num_threads);
int omp_get_dynamic (void);
void omp_set_nested (int nested);
int omp_get_nested (void);
#endif

/* da_string.c */
char *da_strchr_replace (char *const str, const char *const fromlist,
		const char *const tolist);
char *da_strtprune (char *const str, const char *const rmlist);
char *da_strhprune (char *const str, const char *const rmlist);
char *da_strtoupper (char *const str);
char *da_strtolower (char *const str);
char *da_strdup (const char *const orgstr);
int da_strcasecmp (const char *const s1, const char *const s2);
int da_strrcmp (const char *const s1, const char *const s2);
char *da_time2str (time_t time);
char *da_getStringKey (const da_StringMap_t * const strmap, const char id);
int da_getStringID (const da_StringMap_t * const strmap,
		const char *const key);

/* da_random.c */
void da_randinit (uint64_t);
uint64_t da_randint64 (void);
uint32_t da_randint32 (void);

/* memory.c */
void da_AllocMatrix (void ***, const size_t elmlen, const size_t ndim1,
		const size_t ndim2);
void da_FreeMatrix (void ***, const size_t, const size_t);
int da_malloc_init ();
void da_malloc_cleanup (const int showstats);
void *da_malloc (const size_t nbytes, const char *const msg);
void *da_realloc (void *oldptr, const size_t nbytes, const char *const msg);
void da_free (void **ptr1, ...);
size_t da_GetCurMemoryUsed ();
size_t da_GetMaxMemoryUsed ();
void da_GetVMInfo (size_t * vmsize, size_t * vmrss);


/* da_io.c */
FILE *da_fopen (const char *const fname, const char *const mode,
		const char *const msg);
void da_fclose (FILE * stream);
ssize_t da_read (const int fd, void *vbuf, const size_t count);
ssize_t da_getline (char **lineptr, size_t * n, FILE * stream);
char **da_readfile (const char *const fname, size_t * r_nlines);
int32_t *da_i32readfile (const char *const fname, size_t * r_nlines);
int64_t *da_i64readfile (const char *const fname, size_t * r_nlines);
int32_t *da_i32readfilebin (const char *const fname, size_t * r_nelmnts);
size_t da_i32writefilebin (const char *const fname, const size_t n,
		const int32_t * const a);
int64_t *da_i64readfilebin (const char *const fname, size_t * r_nelmnts);
size_t da_i64writefilebin (const char *const fname, const size_t n,
		const int64_t * const a);
float *da_freadfilebin (const char *const fname, size_t * r_nelmnts);
size_t da_fwritefilebin (const char *const fname, const size_t n,
		const float *const a);
double *da_dreadfilebin (const char *const fname, size_t * r_nelmnts);
size_t da_dwritefilebin (const char *const fname, const size_t n,
		const double *const a);

/* da_fs.c */
int da_fexists (const char *const fname);
int da_dexists (const char *const fname);
ssize_t da_getfsize (const char *const fname);
void da_getfilestats (const char *const fname, size_t * r_nlines,
		size_t * r_ntokens, size_t * r_max_nlntokens,
		size_t * r_nbytes);
char *da_getbasename (const char *const path);
char *da_getextname (const char *const path);
char *da_getfilename (const char *const path);
char *da_getpathname (const char *const path);
int da_mkpath (const char *const path);
int da_rmpath (const char *const path);

/* da_mcore.c */
da_mcore_t *da_mcoreCreate (const size_t coresize);
da_mcore_t *da_gkmcoreCreate (void);
void da_mcoreDestroy (da_mcore_t ** r_mcore, const int showstats);
void da_gkmcoreDestroy (da_mcore_t ** r_mcore, const int showstats);
void *da_mcoreMalloc (da_mcore_t * const mcore, size_t nbytes);
void da_mcorePush (da_mcore_t * const mcore);
void da_gkmcorePush (da_mcore_t * const mcore);
void da_mcorePop (da_mcore_t * const mcore);
void da_gkmcorePop (da_mcore_t * const mcore);
void da_mcoreAdd (da_mcore_t * const mcore, const int type,
		const size_t nbytes, void *const ptr);
void da_gkmcoreAdd (da_mcore_t * const mcore, const int type,
		const size_t nbytes, void *const ptr);
void da_mcoreDel (da_mcore_t * const mcore, void *const ptr);
void da_gkmcoreDel (da_mcore_t * const mcore, void *const ptr);

/* csr.c  */
da_csr_t *da_csr_Create ();
void da_csr_Init (da_csr_t * const mat);
void da_csr_Free (da_csr_t ** const mat);
void da_csr_FreeAll (da_csr_t ** ptr1, ...);
void da_csr_FreeBase (da_csr_t * const mat, const char type);
void da_csr_LoadBases (da_csr_t * const csr);
void da_csr_FreeContents (da_csr_t * const mat);
da_csr_t *da_csr_Copy (const da_csr_t * const mat);
da_csr_t * da_csr_Join (const da_csr_t * const mat1, const da_csr_t * const mat2);
da_csr_t *da_csr_ExtractSubmatrix (const da_csr_t * const mat,
		const idx_t rstart, const idx_t nrows);
da_csr_t *da_csr_ExtractRows (const da_csr_t * const mat, const idx_t nrows,
		const idx_t * const rind);
void da_csr_ExtractRowsInto (const da_csr_t * const mat,
		da_csr_t * const nmat, const idx_t nrows,
		const idx_t * const rind);
da_csr_t *da_csr_ExtractPartition (const da_csr_t * const mat,
		const idx_t * const part,
		const idx_t pid);
da_csr_t **da_csr_Split (const da_csr_t * const mat,
		const idx_t * const color);
da_csr_t *da_csr_Read (const char *const filename, const char format,
		char readvals, char numbering);
void da_csr_Write (const da_csr_t * const mat, const char *const filename,
		const char format, char writevals, char numbering);
void da_csr_PrintInfo (const da_csr_t * const mat, const char *const name,
		const char *const suffix);
void da_csr_Print (const da_csr_t * const mat);
char da_csr_isClutoOrCsr (const char *const file);
da_csr_t *da_csr_Prune (const da_csr_t * const mat, const char what,
		const idx_t minf, const idx_t maxf);
da_csr_t *da_csr_LowFilter (const da_csr_t * const mat, const char what,
		const char norm, const float fraction);
da_csr_t *da_csr_HighAvgFilter (const da_csr_t * const mat, const char what,
		const float percent);
da_csr_t *da_csr_topKPlusFilter (const da_csr_t * const mat,
		const char what, const idx_t topk,
		const val_t keepval);
da_csr_t *da_csr_ZScoreFilter (const da_csr_t * const mat, const char what,
		const float zscore);
void da_csr_CompactColumns (da_csr_t * const mat);
void da_csr_CompactRows (da_csr_t * const mat);
void da_csr_SortIndices (da_csr_t * const mat, const char what);
char da_csr_CheckSortedIndex (da_csr_t * const mat, const char what);
void da_csr_SortValues (da_csr_t * const mat, const char what,
		const char how);
void da_csr_CreateIndex (da_csr_t * const mat, const char what);
void da_csr_Normalize (da_csr_t * const mat, const char what,
		const char norm);
void da_csr_Scale (da_csr_t * const mat, const char type);
void da_csr_ComputeSums (da_csr_t * const mat, const char what);
double *da_csr_ComputeMeans (const da_csr_t * const mat, const char what);
void da_csr_ComputeSquaredNorms (da_csr_t * const mat, const char what);
val_t da_csr_ComputeSimilarity (const da_csr_t * const mat,
		const idx_t rc1, const idx_t rc2,
		const char what, const char simtype);
char da_csr_Compare (const da_csr_t * const a, const da_csr_t * const b,
		const double p);
void da_csr_Grow (da_csr_t * const mat, const ptr_t newNnz);
val_t da_csr_partialDotProduct (const ptr_t * const rowptr,
		const ptr_t * const endptr,
		const idx_t * const rowind,
		const val_t * const rowval, const idx_t a,
		const idx_t b);
val_t da_csr_dotProduct (const ptr_t * const rowptr,
		const idx_t * const rowind,
		const val_t * const rowval, const idx_t a,
		const idx_t b);
idx_t da_csr_GetSimilarSmallerRows (const da_csr_t * const mat,
		const idx_t rid, const char noSelfSim,
		const idx_t nqterms,
		const idx_t * const qind,
		const val_t * const qval,
		const char simtype, idx_t nsim,
		const float minsim, da_ivkv_t * hits,
		idx_t * i_marker, da_ivkv_t * i_cand);
void da_csr_Transfer (da_csr_t * from, da_csr_t * to);


/* select.c */
idx_t da_ivkvkselectd (size_t n, idx_t topk, da_ivkv_t * cand);
idx_t da_ivkvkselecti (size_t n, idx_t topk, da_ivkv_t * cand);
idx_t da_iakvkselectd (size_t n, idx_t topk, da_iakv_t * cand);
idx_t da_iakvkselecti (size_t n, idx_t topk, da_iakv_t * cand);


/* sort.c */

#define DA_MKSORT_PROTO(PRFX, TYPE) \
		void     PRFX ## sorti(size_t n, TYPE *base);\
		void     PRFX ## sortd(size_t n, TYPE *base);\

DA_MKSORT_PROTO (da_p, ptr_t)
DA_MKSORT_PROTO (da_i, idx_t)
DA_MKSORT_PROTO (da_v, val_t)
DA_MKSORT_PROTO (da_a, accum_t)
DA_MKSORT_PROTO (da_il, da_il_t *)
DA_MKSORT_PROTO (da_bil, da_bil_t *)
DA_MKSORT_PROTO (da_ppkv, da_ppkv_t)
DA_MKSORT_PROTO (da_pikv, da_pikv_t)
DA_MKSORT_PROTO (da_pvkv, da_pvkv_t)
DA_MKSORT_PROTO (da_pckv, da_pckv_t)
DA_MKSORT_PROTO (da_pi32kv, da_pi32kv_t)
DA_MKSORT_PROTO (da_pi64kv, da_pi64kv_t)
DA_MKSORT_PROTO (da_pzkv, da_pzkv_t)
DA_MKSORT_PROTO (da_pfkv, da_pfkv_t)
DA_MKSORT_PROTO (da_pdkv, da_pdkv_t)
DA_MKSORT_PROTO (da_pskv, da_pskv_t)
DA_MKSORT_PROTO (da_ipkv, da_ipkv_t)
DA_MKSORT_PROTO (da_iikv, da_iikv_t)
DA_MKSORT_PROTO (da_ivkv, da_ivkv_t)
DA_MKSORT_PROTO (da_ickv, da_ickv_t)
DA_MKSORT_PROTO (da_ii32kv, da_ii32kv_t)
DA_MKSORT_PROTO (da_ii64kv, da_ii64kv_t)
DA_MKSORT_PROTO (da_izkv, da_izkv_t)
DA_MKSORT_PROTO (da_ifkv, da_ifkv_t)
DA_MKSORT_PROTO (da_idkv, da_idkv_t)
DA_MKSORT_PROTO (da_iakv, da_iakv_t)
DA_MKSORT_PROTO (da_iskv, da_iskv_t)
DA_MKSORT_PROTO (da_upkv, da_upkv_t)
DA_MKSORT_PROTO (da_uikv, da_uikv_t)
DA_MKSORT_PROTO (da_uvkv, da_uvkv_t)
DA_MKSORT_PROTO (da_uakv, da_uakv_t)
DA_MKSORT_PROTO (da_uckv, da_uckv_t)
DA_MKSORT_PROTO (da_ui32kv, da_ui32kv_t)
DA_MKSORT_PROTO (da_ui64kv, da_ui64kv_t)
DA_MKSORT_PROTO (da_uzkv, da_uzkv_t)
DA_MKSORT_PROTO (da_ufkv, da_ufkv_t)
DA_MKSORT_PROTO (da_udkv, da_udkv_t)
DA_MKSORT_PROTO (da_uskv, da_uskv_t) DA_MKSORT_PROTO (da_uukv, da_uukv_t)
/**
 * Memory allocation functions
 */
 // base types: ptr, idx, and val

 DA_MKALLOC (da_p, ptr_t)
 DA_MKALLOC (da_i, idx_t)
 DA_MKALLOC (da_v, val_t)
 DA_MKALLOC (da_a, accum_t)
 DA_MKALLOC (da_c, char)
 DA_MKALLOC (da_uc, unsigned char)
 DA_MKALLOC (da_i8, int8_t)
 DA_MKALLOC (da_ui8, uint8_t)
 DA_MKALLOC (da_i16, int16_t)
 DA_MKALLOC (da_ui16, uint16_t)
 DA_MKALLOC (da_i32, int32_t)
 DA_MKALLOC (da_ui32, uint32_t)
 DA_MKALLOC (da_i64, int64_t)
 DA_MKALLOC (da_ui64, uint64_t)
 DA_MKALLOC (da_u, uint)
 DA_MKALLOC (da_z, ssize_t)
 DA_MKALLOC (da_uz, size_t)
 DA_MKALLOC (da_f, float)
 DA_MKALLOC (da_d, double)
 DA_MKALLOC (da_l, long)
 DA_MKALLOC (da_ul, unsigned long) DA_MKALLOC (da_s, da_sim_t)
 // ptr_t key, various vals

 DA_MKALLOC (da_ppkv, da_ppkv_t)
 DA_MKALLOC (da_pikv, da_pikv_t)
 DA_MKALLOC (da_pvkv, da_pvkv_t)
 DA_MKALLOC (da_pckv, da_pckv_t)
 DA_MKALLOC (da_pi32kv, da_pi32kv_t)
 DA_MKALLOC (da_pi64kv, da_pi64kv_t)
 DA_MKALLOC (da_pzkv, da_pzkv_t)
 DA_MKALLOC (da_pfkv, da_pfkv_t)
 DA_MKALLOC (da_pdkv, da_pdkv_t) DA_MKALLOC (da_pskv, da_pskv_t)
 // idx_t key, various vals

 DA_MKALLOC (da_ipkv, da_ipkv_t)
 DA_MKALLOC (da_iikv, da_iikv_t)
 DA_MKALLOC (da_ivkv, da_ivkv_t)
 DA_MKALLOC (da_iakv, da_iakv_t)
 DA_MKALLOC (da_ickv, da_ickv_t)
 DA_MKALLOC (da_ii32kv, da_ii32kv_t)
 DA_MKALLOC (da_ii64kv, da_ii64kv_t)
 DA_MKALLOC (da_izkv, da_izkv_t)
 DA_MKALLOC (da_ifkv, da_ifkv_t)
 DA_MKALLOC (da_idkv, da_idkv_t) DA_MKALLOC (da_iskv, da_iskv_t)
 // size_t key, various vals

 DA_MKALLOC (da_upkv, da_upkv_t)
 DA_MKALLOC (da_uikv, da_uikv_t)
 DA_MKALLOC (da_uvkv, da_uvkv_t)
 DA_MKALLOC (da_uakv, da_uakv_t)
 DA_MKALLOC (da_uckv, da_uckv_t)
 DA_MKALLOC (da_ui32kv, da_ui32kv_t)
 DA_MKALLOC (da_ui64kv, da_ui64kv_t)
 DA_MKALLOC (da_uzkv, da_uzkv_t)
 DA_MKALLOC (da_ufkv, da_ufkv_t)
 DA_MKALLOC (da_udkv, da_udkv_t)
 DA_MKALLOC (da_uskv, da_uskv_t) DA_MKALLOC (da_uukv, da_uukv_t)
 // iiakv
 DA_MKALLOC (da_iiakv, da_iiakv_t)
 // l2ap_hash_t
 DA_MKALLOC (da_l2aph, l2ap_hash_t)
 /* pqueue.c */
 da_iapq_t *da_iapqCreate (size_t maxsize, size_t maxnodes);
da_iapq_t *da_iapqCreateShared (size_t maxsize, size_t maxnodes,
		ptr_t * locator);
void da_iapqInit (da_iapq_t * queue, size_t maxsize, size_t maxnodes,
		ptr_t * locator);
void da_iapqInitLocator (da_iapq_t * queue);
void da_iapqResetLocator (da_iapq_t * queue);
void da_iapqReset (da_iapq_t * queue);
void da_iapqFree (da_iapq_t * queue);
void da_iapqDestroy (da_iapq_t * queue);
size_t da_iapqLength (da_iapq_t * queue);
int da_iapqInsert (da_iapq_t * queue, idx_t node, accum_t val);
int da_iapqInsertHeap (da_iapq_t * queue, idx_t node, val_t val);
void da_iapqUpdate (da_iapq_t * queue, idx_t node, accum_t newval);
void da_iapqUpdateHeap (da_iapq_t * queue, ptr_t loc, idx_t node,
		val_t newval);
int da_iapqDelete (da_iapq_t * queue, idx_t node);
da_iakv_t da_iapqGetTop (da_iapq_t * queue);
da_iakv_t da_iapqGetTopHeap (da_iapq_t * queue);
idx_t da_iapqSeeTopKey (da_iapq_t * queue);
accum_t da_iapqSeeTopVal (da_iapq_t * queue);
accum_t da_iapqSeeVal (da_iapq_t * queue, idx_t node);
accum_t da_iapqSeeValHeap (da_iapq_t * queue, idx_t node);
int da_iapqExists (da_iapq_t * queue, idx_t node);
int da_iapqExistsHeap (da_iapq_t * queue, idx_t node);
int da_iapqCheckHeap (da_iapq_t * queue);
void da_iapqPrintAll (da_iapq_t * queue);

DA_MKPQUEUE_PROTO (da_ivpq, da_ivpq_t, idx_t, val_t)
DA_MKPQUEUE_PROTO (da_ivmq, da_ivpq_t, idx_t, val_t)
/**
 * BLAS functions
 */

DA_MKBLAS (da_p, da_up, ptr_t, ptr_t)
DA_MKBLAS (da_i, da_ui, idx_t, idx_t)
DA_MKBLAS (da_v, da_uv, val_t, val_t)
DA_MKBLAS (da_a, da_ua, accum_t, accum_t)
DA_MKBLAS (da_c, da_uc, char, int)
DA_MKBLAS (da_i32, da_ui32, int32_t, int32_t)
DA_MKBLAS (da_i64, da_ui64, int64_t, int64_t)
DA_MKBLAS (da_z, da_uz, ssize_t, ssize_t)
DA_MKBLAS (da_f, da_uf, float, float)
DA_MKBLAS (da_d, da_ud, double, double)
DA_MKBLAS (da_u, da_uu, uint, uint)
/**
  * RAND functions
  */

DA_MKRANDOM (da_p, size_t, ptr_t)
DA_MKRANDOM (da_i, size_t, idx_t)
DA_MKRANDOM (da_v, size_t, val_t)
DA_MKRANDOM (da_a, size_t, accum_t)
DA_MKRANDOM (da_c, size_t, char)
DA_MKRANDOM (da_f, size_t, float)
DA_MKRANDOM (da_d, size_t, double)
DA_MKRANDOM (da_z, size_t, ssize_t) DA_MKRANDOM (da_u, size_t, uint)
#ifdef __cplusplus
}
#endif

#endif
