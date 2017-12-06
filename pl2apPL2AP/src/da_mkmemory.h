/*!
\file  da_mkmemory.h
\brief Templates for memory allocation routines

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Aug 2013.

\date   Started 3/29/07
\author George
*/

#ifndef _COSGRAPH_MKMEMORY_H_
#define _COSGRAPH_MKMEMORY_H_


#define DA_MKALLOC(PRFX, TYPE)\
/*************************************************************************/\
/*! The macro for da_?malloc()-class of routines */\
/**************************************************************************/\
static TYPE* PRFX ## malloc(const size_t n, const char* const msg)\
{\
  return (TYPE *)da_malloc(sizeof(TYPE)*n, msg);\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?realloc()-class of routines */\
/**************************************************************************/\
static TYPE* PRFX ## realloc(TYPE *ptr, const size_t n, const char* const msg)\
{\
  return (TYPE *)da_realloc((void *)ptr, sizeof(TYPE)*n, msg);\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?set()-class of routines */\
/*************************************************************************/\
static TYPE* PRFX ## set(const size_t n, const TYPE val, TYPE* const x)\
{\
  size_t i;\
\
  for (i=0; i<n; i++)\
    x[i] = val;\
\
  return x;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?set()-class of routines */\
/*************************************************************************/\
static TYPE* PRFX ## setzero(const size_t n, TYPE* const x)\
{\
  memset(x, 0, sizeof(TYPE)*n); \
\
  return x;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?set()-class of routines */\
/*************************************************************************/\
static TYPE* PRFX ## copy(const size_t n, TYPE *a, TYPE *b)\
{\
  return (TYPE *)memmove((void *)b, (void *)a, sizeof(TYPE)*n);\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?smalloc()-class of routines */\
/**************************************************************************/\
static TYPE* PRFX ## smalloc(const size_t n, const TYPE ival, const char* msg)\
{\
  TYPE *ptr;\
\
  ptr = (TYPE *)da_malloc(sizeof(TYPE)*n, msg);\
  if (ptr == NULL) \
    return NULL; \
\
  return PRFX ## set(n, ival, ptr); \
}\
\
\
/*************************************************************************/\
/*! The macro for da_?AllocMatrix()-class of routines */\
/**************************************************************************/\
static TYPE** PRFX ## AllocMatrix(const size_t ndim1, const size_t ndim2, const TYPE value, const char* const msg)\
{\
  ssize_t i, j;\
  TYPE **matrix;\
\
  matrix = (TYPE **)da_malloc(ndim1*sizeof(TYPE *), msg);\
  if (matrix == NULL) \
    return NULL;\
\
  for (i=0; i<ndim1; i++) { \
    matrix[i] = PRFX ## smalloc(ndim2, value, msg);\
    if (matrix[i] == NULL) { \
      for (j=0; j<i; j++) \
        da_free((void **)&matrix[j], LTERM); \
      return NULL; \
    } \
  }\
\
  return matrix;\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?AllocMatrix()-class of routines */\
/**************************************************************************/\
static void PRFX ## FreeMatrix(TYPE ***r_matrix, const size_t ndim1, const size_t ndim2)\
{\
  ssize_t i;\
  TYPE **matrix;\
\
  if (*r_matrix == NULL) \
    return; \
\
  matrix = *r_matrix;\
\
  for (i=0; i<ndim1; i++) \
    da_free((void **)&(matrix[i]), LTERM);\
\
  da_free((void **)r_matrix, LTERM);\
}\
\
\
/*************************************************************************/\
/*! The macro for da_?SetMatrix()-class of routines */\
/**************************************************************************/\
static void PRFX ## SetMatrix(TYPE **matrix, const size_t ndim1, const size_t ndim2, const TYPE value)\
{\
  idx_t i, j;\
\
  for (i=0; i<ndim1; i++) {\
    for (j=0; j<ndim2; j++)\
      matrix[i][j] = value;\
  }\
}\



#endif
