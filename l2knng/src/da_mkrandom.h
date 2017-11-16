/*!
\file  da_mkrandom.h
\brief Templates for portable random number generation

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Aug 2013.

\date   Started 5/17/07
\author George
*/


#ifndef _L2KNNG_MKRANDOM_H
#define _L2KNNG_MKRANDOM_H

/*************************************************************************/\
/*! The generator for the rand() related routines.  \
   \params RNGT  the datatype that defines the range of values over which\
                 random numbers will be generated\
   \params VALT  the datatype that defines the contents of the array to \
                 be permuted by randArrayPermute() \
   \params FPRFX the function prefix \
*/\
/**************************************************************************/\
#define DA_MKRANDOM(FPRFX, RNGT, VALT)\
/*************************************************************************/\
/*! Initializes the generator */ \
/**************************************************************************/\
static void FPRFX ## srand(RNGT seed) \
{\
  da_randinit((uint64_t) seed);\
}\
\
\
/*************************************************************************/\
/*! Returns a random number */ \
/**************************************************************************/\
static RNGT FPRFX ## rand() \
{\
  if (sizeof(RNGT) <= sizeof(int32_t)) \
    return (RNGT)da_randint32(); \
  else \
    return (RNGT)da_randint64(); \
}\
\
\
/*************************************************************************/\
/*! Returns a random number between [0, max) */ \
/**************************************************************************/\
static RNGT FPRFX ## randInRange(const RNGT max) \
{\
  return (RNGT)((FPRFX ## rand())%max); \
}\
\
\
/*************************************************************************/\
/*! Randomly permutes the elements of an array p[]. \
    flag == 1, p[i] = i prior to permutation, \
    flag == 0, p[] is not initialized. */\
/**************************************************************************/\
static void FPRFX ## randArrayPermute(const RNGT n, VALT* const p, const RNGT nshuffles, const int flag)\
{\
  RNGT i, u, v;\
  VALT tmp;\
\
  if (flag == 1) {\
    for (i=0; i<n; i++)\
      p[i] = (VALT)i;\
  }\
\
  if (n < 10) {\
    for (i=0; i<n; i++) {\
      v = FPRFX ## randInRange(n);\
      u = FPRFX ## randInRange(n);\
      DA_SWAP(p[v], p[u], tmp);\
    }\
  }\
  else {\
    for (i=0; i<nshuffles; i++) {\
      v = FPRFX ## randInRange(n-3);\
      u = FPRFX ## randInRange(n-3);\
      /*DA_SWAP(p[v+0], p[u+0], tmp);*/\
      /*DA_SWAP(p[v+1], p[u+1], tmp);*/\
      /*DA_SWAP(p[v+2], p[u+2], tmp);*/\
      /*DA_SWAP(p[v+3], p[u+3], tmp);*/\
      DA_SWAP(p[v+0], p[u+2], tmp);\
      DA_SWAP(p[v+1], p[u+3], tmp);\
      DA_SWAP(p[v+2], p[u+0], tmp);\
      DA_SWAP(p[v+3], p[u+1], tmp);\
    }\
  }\
}\
\
\
/*************************************************************************/\
/*! Randomly permutes the elements of an array p[]. \
    flag == 1, p[i] = i prior to permutation, \
    flag == 0, p[] is not initialized. */\
/**************************************************************************/\
static void FPRFX ## randArrayPermuteFine(const RNGT n, VALT* const p, const int flag)\
{\
  RNGT i, v;\
  VALT tmp;\
\
  if (flag == 1) {\
    for (i=0; i<n; i++)\
      p[i] = (VALT)i;\
  }\
\
  for (i=0; i<n; i++) {\
    v = FPRFX ## randInRange(n);\
    DA_SWAP(p[i], p[v], tmp);\
  }\
}\



#endif
