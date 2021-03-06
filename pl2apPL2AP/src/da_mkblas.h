/*!
\file  da_mkblas.h
\brief Templates for BLAS-like routines (modified from GKlib)
Note that key value pairs are reversed in this library vs. GKlib

\date   Started 2/27/2013
\author David C. Anastasiu
*/

#ifndef _COSGRAPH_MKBLAS_H_
#define _COSGRAPH_MKBLAS_H_


#define DA_MKBLAS(PRFX, KVPRFX, TYPE, OUTTYPE) \
/*************************************************************************/\
/*! The macro for da_?incset()-class of routines */\
/*************************************************************************/\
static TYPE *PRFX ## incset(const size_t n, const TYPE baseval, TYPE* const x)\
{\
  size_t i;\
\
  for (i=0; i<n; i++)\
    x[i] = baseval+i;\
\
  return x;\
}\
\
/*************************************************************************/\
/*! The macro for da_?max()-class of routines */\
/*************************************************************************/\
static TYPE PRFX ## max(const size_t n, const TYPE* const x)\
{\
  size_t i, max=0; \
\
  if (n <= 0) return (TYPE) 0;\
\
  for (i=1; i<n; i++)\
    max = (x[i] > x[max] ? i : max);\
\
  return x[max];\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?min()-class of routines */\
/*************************************************************************/\
static TYPE PRFX ## min(const size_t n, const TYPE* const x)\
{\
  size_t i, min=0;\
\
  if (n <= 0) return (TYPE) 0;\
\
  for (i=1; i<n; i++)\
    min = (x[i] < x[min] ? i : min);\
\
  return x[min];\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?argmax()-class of routines */\
/*************************************************************************/\
static size_t PRFX ## argmax(const size_t n, const TYPE* const x)\
{\
  size_t i, max=0;\
\
  for (i=1; i<n; i++)\
    max = (x[i] > x[max] ? i : max);\
\
  return max;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?argmin()-class of routines */\
/*************************************************************************/\
static size_t PRFX ## argmin(const size_t n, const TYPE* const x)\
{\
  size_t i, min=0;\
\
  for (i=1; i<n; i++)\
    min = (x[i] < x[min] ? i : min);\
\
  return min;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?sum()-class of routines */\
/**************************************************************************/\
static OUTTYPE PRFX ## sum(const size_t n, const TYPE* x, const size_t incx)\
{\
  size_t i;\
  OUTTYPE sum = 0;\
\
  for (i=0; i<n; i++, x+=incx)\
    sum += (*x);\
\
  return sum;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?count()-class of routines */\
/**************************************************************************/\
static size_t PRFX ## count(const size_t n, const TYPE* const x, const TYPE search)\
{\
  size_t i, cnt=0;\
\
  for (i=0; i<n; i++)\
    if(x[i] == search) cnt++;\
\
  return cnt;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?add()-class of routines - add to each element */\
/**************************************************************************/\
static TYPE *PRFX ## add(const size_t n, const TYPE alpha, TYPE* x, const size_t incx)\
{\
  size_t i;\
\
  for (i=0; i<n; i++, x+=incx)\
    (*x) += alpha;\
\
  return x;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?scale()-class of routines */\
/**************************************************************************/\
static TYPE *PRFX ## scale(const size_t n, const TYPE alpha, TYPE* x, const size_t incx)\
{\
  size_t i;\
\
  for (i=0; i<n; i++, x+=incx)\
    (*x) *= alpha;\
\
  return x;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?norm2()-class of routines */\
/**************************************************************************/\
static OUTTYPE PRFX ## norm2(const size_t n, const TYPE* x, const size_t incx)\
{\
  size_t i;\
  OUTTYPE partial = 0;\
\
  for (i=0; i<n; i++, x+=incx)\
    partial += (*x) * (*x);\
\
  return (partial > 0 ? (OUTTYPE)sqrt((double)partial) : (OUTTYPE)0);\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?dot()-class of routines */\
/**************************************************************************/\
static OUTTYPE PRFX ## dot(const size_t n, const TYPE* x, const size_t incx, const TYPE* y, const size_t incy)\
{\
  size_t i;\
  OUTTYPE partial = 0.0;\
 \
  for (i=0; i<n; i++, x+=incx, y+=incy)\
    partial += (*x) * (*y);\
\
  return partial;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?arreq()-class of routines */\
/**************************************************************************/\
static char PRFX ## arreq(const size_t n, const TYPE* const x, const TYPE* const y)\
{\
  size_t i;\
 \
  for (i=0; i<n; i++)\
    if(x[i] != y[i]) return 0;\
\
  return 1;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?arreq()-class of routines */\
/**************************************************************************/\
static size_t PRFX ## arreqcnt(const size_t n, const TYPE* const x, const TYPE* const y)\
{\
  size_t i, j;\
 \
  for (j=0, i=0; i<n; i++)\
    if(x[i] != y[i]) j++;\
\
  return j;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?arreq_p()-class of routines */\
/**************************************************************************/\
static char PRFX ## arreq_p(const size_t n, const TYPE* const x, const TYPE* const y, const double p)\
{\
  size_t i;\
 \
  for (i=0; i<n; i++)\
    if(da_abs(x[i] - y[i]) > p) return 0;\
\
  return 1;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?axpy()-class of routines */\
/**************************************************************************/\
static TYPE *PRFX ## axpy(const size_t n, const TYPE alpha, const TYPE* x, const size_t incx, TYPE* y, const size_t incy)\
{\
  size_t i;\
  TYPE *y_in = y;\
\
  for (i=0; i<n; i++, x+=incx, y+=incy)\
    *y += alpha*(*x);\
\
  return y_in;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?argmax_n()-class of routines */\
/*************************************************************************/\
static size_t PRFX ## argmax_n(const size_t n, const TYPE* const x, const size_t k)\
{\
  size_t i, max_n;\
  KVPRFX ## kv_t *cand;\
\
  cand = KVPRFX ## kvmalloc(n, "DA_ARGMAX_N: cand");\
\
  for (i=0; i<n; i++) {\
    cand[i].key = i;\
    cand[i].val = x[i];\
  }\
  KVPRFX ## kvsortd(n, cand);\
\
  max_n = cand[k-1].key;\
\
  da_free((void *)&cand, LTERM);\
\
  return max_n;\
}\


#define DA_MKBLAS_PROTO(PRFX, TYPE, OUTTYPE) \
  TYPE    *PRFX ## incset(size_t n, TYPE baseval, TYPE *x);\
  TYPE     PRFX ## max(size_t n, TYPE *x);\
  TYPE     PRFX ## min(size_t n, TYPE *x);\
  size_t   PRFX ## argmax(size_t n, TYPE *x);\
  size_t   PRFX ## argmin(size_t n, TYPE *x);\
  OUTTYPE  PRFX ## sum(size_t n, TYPE *x, size_t incx);\
  TYPE    *PRFX ## scale(size_t n, TYPE alpha, TYPE *x, size_t incx);\
  OUTTYPE  PRFX ## norm2(size_t n, TYPE *x, size_t incx);\
  OUTTYPE  PRFX ## dot(size_t n, TYPE *x, size_t incx, TYPE *y, size_t incy);\
  char     PRFX ## arreq(size_t n, TYPE *x, TYPE *y);\
  char     PRFX ## arreq_p(size_t n, TYPE *x, TYPE *y, double p);\
  TYPE    *PRFX ## axpy(size_t n, TYPE alpha, TYPE *x, size_t incx, TYPE *y, size_t incy);\
  size_t   PRFX ## argmax_n(size_t n, TYPE *x, size_t k);\



#endif
