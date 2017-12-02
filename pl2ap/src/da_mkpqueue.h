/*!
\file  da_mkpqueue.h
\brief Templates for priority queues

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Oct 2014.

\date   Started 4/09/07
\author George
\version\verbatim $Id: da_mkpqueue.h 16221 2014-02-15 20:33:26Z karypis $ \endverbatim
*/


#ifndef _COSGRAPH_MKPQUEUE_H
#define _COSGRAPH_MKPQUEUE_H


#define DA_MKPQUEUE(FPRFX, PQT, KVT, KT, VT, KVMALLOC, VMAX, VAL_LT)\
/****************************************************************************/\
/*! This function creates and initializes a size-constrained priority queue */\
/****************************************************************************/\
PQT *FPRFX ## Create(size_t maxsize, size_t maxnodes)\
{\
  PQT *queue; \
\
  queue = (PQT *)da_malloc(sizeof(PQT), "da_pqCreate: queue");\
  FPRFX ## Init(queue, maxsize, maxnodes, NULL);\
\
  return queue;\
}\
\
\
/**************************************************************************/\
/*! This function creates and initializes a priority queue using a shared */\
/*! locator array, which should be of size maxnodes and pre-initialized   */\
/*! with -1 values.                                                       */\
/**************************************************************************/\
PQT *FPRFX ## CreateShared(size_t maxsize, size_t maxnodes, ptr_t *locator)\
{\
  PQT *queue; \
\
  queue = (PQT *)da_malloc(sizeof(PQT), "da_pqCreate: queue");\
  FPRFX ## Init(queue, maxsize, maxnodes, locator);\
\
  return queue;\
}\
\
\
/**************************************************************************/\
/*! This function initializes the data structures of the priority queue   */\
/**************************************************************************/\
void FPRFX ## Init(PQT *queue, size_t maxsize, size_t maxnodes, ptr_t *locator)\
{\
  queue->nnodes = 0;\
  queue->maxsize = maxsize;\
  queue->maxnodes = maxnodes;\
\
  queue->heap    = KVMALLOC(maxsize, "da_PQInit: heap");\
  if(locator)\
    queue->locator = locator;\
  else\
    queue->locator = da_psmalloc(maxnodes, -1, "da_PQInit: locator");\
}\
\
\
/**************************************************************************/\
/*! Initialize locator array with values from the heap                    */\
/**************************************************************************/\
void FPRFX ## InitLocator(PQT *queue)\
{\
  ptr_t i;\
  size_t nnodes;\
  ptr_t *locator;\
  KVT *heap;\
\
  heap    = queue->heap;\
  locator = queue->locator;\
  nnodes  = queue->nnodes;\
\
  if (!heap)\
    return;\
\
  for (i=1; i<nnodes; ++i)\
    locator[heap[i].key] = i;\
\
}\
\
\
/**************************************************************************/\
/*! Reset locator array based on values from the heap                     */\
/**************************************************************************/\
void FPRFX ## ResetLocator(PQT *queue)\
{\
  ptr_t i;\
  size_t nnodes;\
  ptr_t *locator;\
  KVT *heap;\
\
  heap    = queue->heap;\
  locator = queue->locator;\
  nnodes  = queue->nnodes;\
\
  if (!heap)\
    return;\
\
  for (i=1; i<nnodes; ++i)\
    locator[heap[i].key] = -1;\
\
}\
\
\
/**************************************************************************/\
/*! This function resets the priority queue                               */\
/**************************************************************************/\
void FPRFX ## Reset(PQT *queue)\
{\
  ptr_t i;\
  ptr_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  for (i=queue->nnodes-1; i>=0; i--)\
    locator[heap[i].key] = -1;\
  queue->nnodes = 0;\
}\
\
\
/***************************************************************************/\
/*! This function frees the internal data structures of the priority queue */\
/***************************************************************************/\
void FPRFX ## Free(PQT *queue)\
{\
  if (queue == NULL) return;\
  da_free((void **)&queue->heap, &queue->locator, LTERM);\
  queue->maxnodes = 0;\
}\
\
\
/*************************************************************************/\
/*! This function frees the internal data structures of the priority queue \
    and the queue itself */\
/**************************************************************************/\
void FPRFX ## Destroy(PQT *queue)\
{\
  if (queue == NULL) return;\
  FPRFX ## Free(queue);\
  da_free((void **)&queue, LTERM);\
}\
\
\
/**************************************************************************/\
/*! This function returns the length of the queue                         */\
/**************************************************************************/\
size_t FPRFX ## Length(PQT *queue)\
{\
  return queue->nnodes;\
}\
\
\
/***************************************************************************/\
/*! This function adds an item in the priority queue if queue is not full. */\
/*! If top item has higher/lower value than the inserted item, it is       */\
/*! replaced by the current item.                                          */\
/***************************************************************************/\
int FPRFX ## Insert(PQT *queue, KT node, VT val)\
{\
  ptr_t i, j;\
  ptr_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  ASSERT(FPRFX ## CheckHeap(queue));\
\
  ASSERT(locator[node] == -1);\
\
/*  if(queue->maxsize == 0) \
    return 0; \
  else if(queue->nnodes == queue->maxsize){  \
    if(VAL_LT(val, heap[0].val)){\
        FPRFX ## GetTop(queue);\
    } else \
      return 0;\
  } */ \
\
  i = queue->nnodes++;\
  while (i > 0) {\
    j = (i-1)>>1;\
    if (VAL_LT(val, heap[j].val)) {\
      heap[i] = heap[j];\
      locator[heap[i].key] = i;\
      i = j;\
    }\
    else\
      break;\
  }\
  ASSERT(i >= 0);\
  heap[i].val   = val;\
  heap[i].key   = node;\
  locator[node] = i;\
\
  ASSERT(FPRFX ## CheckHeap(queue));\
\
  return 1;\
}\
\
\
/*************************************************************************/\
/*! This function deletes an item from the priority queue */\
/**************************************************************************/\
int FPRFX ## Delete(PQT *queue, KT node)\
{\
  ptr_t i, j;\
  size_t nnodes;\
  VT newval, oldval;\
  ptr_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  ASSERT(locator[node] != -1);\
  ASSERT(heap[locator[node]].key == node);\
\
  ASSERT(FPRFX ## CheckHeap(queue));\
\
  i = locator[node];\
  locator[node] = -1;\
\
  if (--queue->nnodes > 0 && heap[queue->nnodes].key != node) {\
    node   = heap[queue->nnodes].key;\
    newval = heap[queue->nnodes].val;\
    oldval = heap[i].val;\
\
    if (VAL_LT(newval, oldval)) { /* Filter-up */\
      while (i > 0) {\
        j = (i-1)>>1;\
        if (VAL_LT(newval, heap[j].val)) {\
          heap[i] = heap[j];\
          locator[heap[i].key] = i;\
          i = j;\
        }\
        else\
          break;\
      }\
    }\
    else { /* Filter down */\
      nnodes = queue->nnodes;\
      while ((j=(i<<1)+1) < nnodes) {\
        if (VAL_LT(heap[j].val, newval)) {\
          if (j+1 < nnodes && VAL_LT(heap[j+1].val, heap[j].val))\
            j++;\
          heap[i] = heap[j];\
          locator[heap[i].key] = i;\
          i = j;\
        }\
        else if (j+1 < nnodes && VAL_LT(heap[j+1].val, newval)) {\
          j++;\
          heap[i] = heap[j];\
          locator[heap[i].key] = i;\
          i = j;\
        }\
        else\
          break;\
      }\
    }\
\
    heap[i].val   = newval;\
    heap[i].key   = node;\
    locator[node] = i;\
  }\
\
  ASSERT(FPRFX ## CheckHeap(queue));\
\
  return 0;\
}\
\
\
/**************************************************************************/\
/*! This function updates the value associated with a particular item     */\
/**************************************************************************/\
void FPRFX ## Update(PQT *queue, KT node, VT newval)\
{\
  ptr_t i, j;\
  size_t nnodes;\
  VT oldval;\
  ptr_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  oldval = heap[locator[node]].val;\
\
  ASSERT(locator[node] != -1);\
  ASSERT(heap[locator[node]].key == node);\
  ASSERT(FPRFX ## CheckHeap(queue));\
\
  i = locator[node];\
\
  if (VAL_LT(newval, oldval)) { /* Filter-up */\
    while (i > 0) {\
      j = (i-1)>>1;\
      if (VAL_LT(newval, heap[j].val)) {\
        heap[i] = heap[j];\
        locator[heap[i].key] = i;\
        i = j;\
      }\
      else\
        break;\
    }\
  }\
  else { /* Filter down */\
    nnodes = queue->nnodes;\
    while ((j=(i<<1)+1) < nnodes) {\
      if (VAL_LT(heap[j].val, newval)) {\
        if (j+1 < nnodes && VAL_LT(heap[j+1].val, heap[j].val))\
          j++;\
        heap[i] = heap[j];\
        locator[heap[i].key] = i;\
        i = j;\
      }\
      else if (j+1 < nnodes && VAL_LT(heap[j+1].val, newval)) {\
        j++;\
        heap[i] = heap[j];\
        locator[heap[i].key] = i;\
        i = j;\
      }\
      else\
        break;\
    }\
  }\
\
  heap[i].val   = newval;\
  heap[i].key   = node;\
  locator[node] = i;\
\
  ASSERT(FPRFX ## CheckHeap(queue));\
\
  return;\
}\
\
\
/*************************************************************************/\
/*! This function returns the item at the top of the queue and removes\
    it from the priority queue */\
/**************************************************************************/\
KT FPRFX ## GetTop(PQT *queue)\
{\
  ptr_t i, j;\
  ptr_t *locator;\
  KVT *heap;\
  KT vtx, node;\
  VT val;\
\
  ASSERT(FPRFX ## CheckHeap(queue));\
\
  if (queue->nnodes == 0)\
    return -1;\
\
  queue->nnodes--;\
\
  heap    = queue->heap;\
  locator = queue->locator;\
\
  vtx = heap[0].key;\
  locator[vtx] = -1;\
\
  if ((i = queue->nnodes) > 0) {\
    val  = heap[i].val;\
    node = heap[i].key;\
    i = 0;\
    while ((j=2*i+1) < queue->nnodes) {\
      if (VAL_LT(heap[j].val, val)) {\
        if (j+1 < queue->nnodes && VAL_LT(heap[j+1].val, heap[j].val))\
          j = j+1;\
        heap[i] = heap[j];\
        locator[heap[i].key] = i;\
        i = j;\
      }\
      else if (j+1 < queue->nnodes && VAL_LT(heap[j+1].val, val)) {\
        j = j+1;\
        heap[i] = heap[j];\
        locator[heap[i].key] = i;\
        i = j;\
      }\
      else\
        break;\
    }\
\
    heap[i].val   = val;\
    heap[i].key   = node;\
    locator[node] = i;\
  }\
\
  ASSERT(FPRFX ## CheckHeap(queue));\
  return vtx;\
}\
\
\
/*************************************************************************/\
/*! This function returns the item at the top of the queue. The item is not\
    deleted from the queue. */\
/**************************************************************************/\
KT FPRFX ## SeeTopKey(PQT *queue)\
{\
  return (queue->nnodes == 0 ? -1 : queue->heap[0].key);\
}\
\
\
/*************************************************************************/\
/*! This function returns the value of the top item. The item is not\
    deleted from the queue. */\
/**************************************************************************/\
VT FPRFX ## SeeTopVal(PQT *queue)\
{\
  return (queue->nnodes == 0 ? VMAX : queue->heap[0].val);\
}\
\
\
/*************************************************************************/\
/*! This function returns the value of a specific item */\
/**************************************************************************/\
VT FPRFX ## SeeVal(PQT *queue, KT node)\
{\
  ptr_t *locator;\
  KVT *heap;\
\
  heap    = queue->heap;\
  locator = queue->locator;\
\
  return heap[locator[node]].val;\
}\
\
\
/*************************************************************************/\
/*! This function returns the value of a specific item */\
/**************************************************************************/\
int FPRFX ## Exists(PQT *queue, KT node)\
{\
  ptr_t *locator;\
\
  locator = queue->locator;\
\
  return locator[node] != -1;\
}\
\
\
/*************************************************************************/\
/*! This functions checks the consistency of the heap */\
/**************************************************************************/\
int FPRFX ## CheckHeap(PQT *queue)\
{\
  ptr_t i, j;\
  size_t nnodes;\
  ptr_t *locator;\
  KVT *heap;\
\
  heap    = queue->heap;\
  locator = queue->locator;\
  nnodes  = queue->nnodes;\
\
  if (!heap)\
    da_errexit("CheckHeap: Heap not initialized!!!");\
\
  if (nnodes == 0)\
    return 1;\
\
  ASSERT(locator[heap[0].key] == 0);\
  for (i=1; i<nnodes; ++i) {\
    ASSERT(locator[heap[i].key] == i);\
    ASSERT(!VAL_LT(heap[i].val, heap[(i-1)/2].val));\
    ASSERT(!VAL_LT(heap[i].val, heap[0].val));\
  }\
\
  for (j=i=0; i<queue->maxnodes; ++i) {\
    if (locator[i] != -1)\
      j++;\
  }\
  ASSERTP(j == nnodes, ("%jd %jd\n", (intmax_t)j, (intmax_t)nnodes));\
\
  return 1;\
}\


#define DA_MKPQUEUE_PROTO(FPRFX, PQT, KT, VT)\
  PQT *  FPRFX ## Create(size_t maxsize, size_t maxnodes);\
  PQT *  FPRFX ## CreateShared(size_t maxsize, size_t maxnodes, ptr_t *locator);\
  void   FPRFX ## Init(PQT *queue, size_t maxsize, size_t maxnodes, ptr_t *locator);\
  void   FPRFX ## InitLocator(PQT *queue);\
  void   FPRFX ## ResetLocator(PQT *queue);\
  void   FPRFX ## Reset(PQT *queue);\
  void   FPRFX ## Free(PQT *queue);\
  void   FPRFX ## Destroy(PQT *queue);\
  size_t FPRFX ## Length(PQT *queue);\
  int    FPRFX ## Insert(PQT *queue, KT node, VT val);\
  int    FPRFX ## Delete(PQT *queue, KT node);\
  void   FPRFX ## Update(PQT *queue, KT node, VT newval);\
  KT     FPRFX ## GetTop(PQT *queue);\
  KT     FPRFX ## SeeTopKey(PQT *queue);\
  VT     FPRFX ## SeeTopVal(PQT *queue);\
  VT     FPRFX ## SeeVal(PQT *queue, KT node);\
  int    FPRFX ## Exists(PQT *queue, KT node);\
  int    FPRFX ## CheckHeap(PQT *queue);\
  void   FPRFX ## PrintAll(PQT *queue);\

/***
 * This is how these macros are used
 *
 * in proto.h:
 * DA_MKPQUEUE_PROTO(da_ivpq, da_ivpq_t, idx_t, val_t)
 *
 * in pqueue.c:
 * #define val_lt(a, b) ((a) < (b))
 * DA_MKPQUEUE(da_ivpq, da_ivpq_t, da_ivkv_t, idx_t, val_t, da_ivkvmalloc, FLT_MAX, val_lt)
 * #undef val_lt
 ***/


#endif
