/*!
\file  pqueue.c
\brief This file implements various max-priority queues.

The priority queues are generated using the GK_MKPQUEUE macro.

\date   Started 3/27/2007
\author George
\version\verbatim $Id: pqueue.c 10711 2011-08-31 22:23:04Z karypis $ \endverbatim
 */

#include "includes.h"


/*************************************************************************/
/*! Create the various priority queues */
/*************************************************************************/


#define val_gt(a, b) ((a) > (b))

/** NOTE: Inserts are not checked in da_ivpq. Check item does not exist before inserting */
DA_MKPQUEUE(da_ivmq, da_ivpq_t, da_ivkv_t, idx_t, val_t, da_ivkvmalloc, -FLT_MAX, val_gt)

#undef val_gt



#define val_lt(a, b) ((a) < (b))

/** NOTE: Inserts are not checked in da_ivpq. Check item does not exist before inserting */
DA_MKPQUEUE(da_ivpq, da_ivpq_t, da_ivkv_t, idx_t, val_t, da_ivkvmalloc, FLT_MAX, val_lt)

/****************************************************************************/
/*! This function creates and initializes a size-constrained priority queue */
/*  in the range [0, inf)                                                   */
/****************************************************************************/
da_iapq_t *da_iapqCreate(size_t maxsize, size_t maxnodes)
{
    da_iapq_t *queue;

    queue = (da_iapq_t *)da_malloc(sizeof(da_iapq_t), "da_pqCreate: queue");
    da_iapqInit(queue, maxsize, maxnodes, NULL);

    return queue;
}


/**************************************************************************/
/*! This function creates and initializes a priority queue using a shared */
/*! locator array, which should be of size maxnodes and pre-initialized   */
/*! with -1 values.                                                       */
/**************************************************************************/
da_iapq_t *da_iapqCreateShared(size_t maxsize, size_t maxnodes, ptr_t *locator)
{
    da_iapq_t *queue;

    queue = (da_iapq_t *)da_malloc(sizeof(da_iapq_t), "da_pqCreate: queue");
    da_iapqInit(queue, maxsize, maxnodes, locator);

    return queue;
}


/**************************************************************************/
/*! This function initializes the data structures of the priority queue   */
/**************************************************************************/
void da_iapqInit(da_iapq_t *queue, size_t maxsize, size_t maxnodes, ptr_t *locator)
{
    queue->nnodes = 0;
    queue->maxsize = maxsize;
    queue->maxnodes = maxnodes;
    ASSERT(maxsize > 0);

    queue->heap    = da_iakvsmalloc(maxsize, (da_iakv_t){-1, 0.0},  "da_PQInit: heap");
    if(locator)
        queue->locator = locator;
    else
        queue->locator = da_psmalloc(maxnodes, -1, "da_PQInit: locator");
}


/**************************************************************************/
/*! Initialize locator array with values from the heap                    */
/**************************************************************************/
void da_iapqInitLocator(da_iapq_t *queue)
{
    ptr_t i;
    size_t nnodes;
    ptr_t *locator;
    da_iakv_t *heap;

    heap    = queue->heap;
    locator = queue->locator;
    nnodes  = queue->nnodes;

    if (!heap)
        return;

    for (i=0; i<nnodes; ++i)
        locator[heap[i].key] = i;

}


/**************************************************************************/
/*! Reset locator array based on values from the heap                     */
/**************************************************************************/
void da_iapqResetLocator(da_iapq_t *queue)
{
    ptr_t i;
    size_t nnodes;
    ptr_t *locator;
    da_iakv_t *heap;

    heap    = queue->heap;
    locator = queue->locator;
    nnodes  = queue->nnodes;

    if (!heap)
        return;

    for (i=0; i<nnodes; ++i)
        locator[heap[i].key] = -1;

}


/**************************************************************************/
/*! This function resets the priority queue                               */
/**************************************************************************/
void da_iapqReset(da_iapq_t *queue)
{
    ptr_t i;
    ptr_t *locator=queue->locator;
    da_iakv_t *heap=queue->heap;

    for (i=queue->nnodes-1; i>=0; i--){
        locator[heap[i].key] = -1;
        heap[i].key = -1;
    }
    queue->nnodes = 0;
}


/***************************************************************************/
/*! This function frees the internal data structures of the priority queue */
/***************************************************************************/
void da_iapqFree(da_iapq_t *queue)
{
    if (queue == NULL) return;
    da_free((void **)&queue->heap, &queue->locator, LTERM);
    queue->maxnodes = 0;
}


/*************************************************************************/
/*! This function frees the internal data structures of the priority queue
    and the queue itself */
/**************************************************************************/
void da_iapqDestroy(da_iapq_t *queue)
{
    if (queue == NULL) return;
    da_iapqFree(queue);
    da_free((void **)&queue, LTERM);
}


/**************************************************************************/
/*! This function returns the length of the queue                         */
/**************************************************************************/
size_t da_iapqLength(da_iapq_t *queue)
{
    return queue->nnodes;
}


/***************************************************************************/
/*! This function adds an item in the priority queue if queue is not full. */
/*! If top item has higher/lower value than the inserted item, it is       */
/*! replaced by the current item.                                          */
/***************************************************************************/
int da_iapqInsert(da_iapq_t *queue, idx_t node, accum_t val)
{
    ptr_t i, j;
    ptr_t *locator=queue->locator;
    da_iakv_t *heap=queue->heap;


    if(locator[node] != -1){
        da_iapqUpdate(queue, node, val);
        return 0;
    }

    if(queue->nnodes == queue->maxsize){
        if(!val_lt(val, heap[0].val)){
            da_iapqGetTop(queue);
        } else
            return 0;
    }

    ASSERT(da_iapqCheckHeap(queue));

    i = queue->nnodes++;
    while (i > 0) {
        j = (i-1)>>1;
        if (val_lt(val, heap[j].val)) {
            heap[i] = heap[j];
            locator[heap[i].key] = i;
            i = j;
        }
        else
            break;
    }
    ASSERT(i >= 0);
    heap[i].val   = val;
    heap[i].key   = node;
    locator[node] = i;

    ASSERT(da_iapqCheckHeap(queue));

    return 1;
}


/***************************************************************************/
/*! This function adds an item in the priority queue if queue is not full. */
/*! If top item has higher/lower value than the inserted item, it is       */
/*! replaced by the current item.                                          */
/*! This version of the function does not make use of the locator array.   */
/***************************************************************************/
int da_iapqInsertHeap(da_iapq_t *queue, idx_t node, val_t val)
{
    ptr_t i, j, loc;
    da_iakv_t *heap=queue->heap;

    /* linear search to find node */
    for(loc = -1, i=0; i < queue->nnodes; ++i)
        if(heap[i].key == node){
            loc = i;
            break;
        }

    if(loc != -1){
        da_iapqUpdateHeap(queue, loc, node, val);
        return 0;
    }

    if(queue->nnodes == queue->maxsize){
        if(!val_lt(val, heap[0].val)){
            da_iapqGetTopHeap(queue);
        } else
            return 0;
    }

    i = queue->nnodes++;
    while (i > 0) {
        j = (i-1)>>1;
        if (val_lt(val, heap[j].val)) {
            heap[i] = heap[j];
            i = j;
        }
        else
            break;
    }
    ASSERT(i >= 0);
    heap[i].val   = val;
    heap[i].key   = node;

    return 1;
}


/**************************************************************************/
/*! This function updates the value associated with a particular item     */
/**************************************************************************/
void da_iapqUpdate(da_iapq_t *queue, idx_t node, accum_t newval)
{
    ptr_t i, j;
    size_t nnodes;
    val_t oldval;
    ptr_t *locator=queue->locator;
    da_iakv_t *heap=queue->heap;

    if(locator[node] == -1)
        return;

    oldval = heap[locator[node]].val;
    if(oldval == newval)
        return;

    ASSERT(heap[locator[node]].key == node);
    ASSERT(da_iapqCheckHeap(queue));

    i = locator[node];

    if (val_lt(newval, oldval)) { /* Filter-up */
        while (i > 0) {
            j = (i-1)>>1;
            if (val_lt(newval, heap[j].val)) {
                heap[i] = heap[j];
                locator[heap[i].key] = i;
                i = j;
            }
            else
                break;
        }
    }
    else { /* Filter down */
        nnodes = queue->nnodes;
        while ((j=(i<<1)+1) < nnodes) {
            if (val_lt(heap[j].val, newval)) {
                if (j+1 < nnodes && val_lt(heap[j+1].val, heap[j].val))
                    j++;
                heap[i] = heap[j];
                locator[heap[i].key] = i;
                i = j;
            }
            else if (j+1 < nnodes && val_lt(heap[j+1].val, newval)) {
                j++;
                heap[i] = heap[j];
                locator[heap[i].key] = i;
                i = j;
            }
            else
                break;
        }
    }

    heap[i].val   = newval;
    heap[i].key   = node;
    locator[node] = i;

    ASSERT(da_iapqCheckHeap(queue));

    return;
}


/**************************************************************************/
/*! This function updates the value associated with a particular item     */
/*! This version of the function does not make use of the locator array.  */
/**************************************************************************/
void da_iapqUpdateHeap(da_iapq_t *queue, ptr_t loc, idx_t node, val_t newval)
{
    ptr_t i, j;
    size_t nnodes;
    val_t oldval;
    da_iakv_t *heap=queue->heap;

    ASSERT(loc != -1)
    oldval = heap[loc].val;
    if(oldval == newval)
        return;
    i = loc;

    if (val_lt(newval, oldval)) { /* Filter-up */
        while (i > 0) {
            j = (i-1)>>1;
            if (val_lt(newval, heap[j].val)) {
                heap[i] = heap[j];
                i = j;
            }
            else
                break;
        }
    }
    else { /* Filter down */
        nnodes = queue->nnodes;
        while ((j=(i<<1)+1) < nnodes) {
            if (val_lt(heap[j].val, newval)) {
                if (j+1 < nnodes && val_lt(heap[j+1].val, heap[j].val))
                    j++;
                heap[i] = heap[j];
                i = j;
            }
            else if (j+1 < nnodes && val_lt(heap[j+1].val, newval)) {
                j++;
                heap[i] = heap[j];
                i = j;
            }
            else
                break;
        }
    }

    heap[i].val   = newval;
    heap[i].key   = node;

    return;
}


/*************************************************************************/
/*! This function deletes an item from the priority queue */
/**************************************************************************/
int da_iapqDelete(da_iapq_t *queue, idx_t node)
{
    ptr_t i, j;
    size_t nnodes;
    val_t newval, oldval;
    ptr_t *locator=queue->locator;
    da_iakv_t *heap=queue->heap;

    ASSERT(locator[node] != -1);
    ASSERT(heap[locator[node]].key == node);

    ASSERT(da_iapqCheckHeap(queue));

    i = locator[node];
    locator[node] = -1;

    if (--queue->nnodes > 0 && heap[queue->nnodes].key != node) {
        node   = heap[queue->nnodes].key;
        newval = heap[queue->nnodes].val;
        oldval = heap[i].val;

        if (val_lt(newval, oldval)) { /* Filter-up */
            while (i > 0) {
                j = (i-1)>>1;
                if (val_lt(newval, heap[j].val)) {
                    heap[i] = heap[j];
                    locator[heap[i].key] = i;
                    i = j;
                }
                else
                    break;
            }
        }
        else { /* Filter down */
            nnodes = queue->nnodes;
            while ((j=(i<<1)+1) < nnodes) {
                if (val_lt(heap[j].val, newval)) {
                    if (j+1 < nnodes && val_lt(heap[j+1].val, heap[j].val))
                        j++;
                    heap[i] = heap[j];
                    locator[heap[i].key] = i;
                    i = j;
                }
                else if (j+1 < nnodes && val_lt(heap[j+1].val, newval)) {
                    j++;
                    heap[i] = heap[j];
                    locator[heap[i].key] = i;
                    i = j;
                }
                else
                    break;
            }
        }

        heap[i].val   = newval;
        heap[i].key   = node;
        locator[node] = i;
    }

    ASSERT(da_iapqCheckHeap(queue));

    return 0;
}


/************************************************************************/
/*! This function returns the item at the top of the queue and removes  */
/*  it from the priority queue                                          */
/************************************************************************/
da_iakv_t da_iapqGetTop(da_iapq_t *queue)
{
    ptr_t i, j;
    ptr_t *locator;
    da_iakv_t *heap;
    idx_t node;
    val_t val;
    da_iakv_t vtx;

    ASSERT(da_iapqCheckHeap(queue));

    if (queue->nnodes == 0)
        return (da_iakv_t){-1,0.0};

    queue->nnodes--;

    heap    = queue->heap;
    locator = queue->locator;

    vtx.key = heap[0].key;
    vtx.val = heap[0].val;
    locator[vtx.key] = -1;

    if ((i = queue->nnodes) > 0) {
        val  = heap[i].val;
        node = heap[i].key;
        i = 0;
        while ((j=2*i+1) < queue->nnodes) {
            if (val_lt(heap[j].val, val)) {
                if (j+1 < queue->nnodes && val_lt(heap[j+1].val, heap[j].val))
                    j = j+1;
                heap[i] = heap[j];
                locator[heap[i].key] = i;
                i = j;
            }
            else if (j+1 < queue->nnodes && val_lt(heap[j+1].val, val)) {
                j = j+1;
                heap[i] = heap[j];
                locator[heap[i].key] = i;
                i = j;
            }
            else
                break;
        }

        heap[i].val   = val;
        heap[i].key   = node;
        locator[node] = i;
    }

    ASSERT(da_iapqCheckHeap(queue));
    return vtx;
}


/*************************************************************************/
/*! This function returns the item at the top of the queue and removes   */
/*  it from the priority queue                                           */
/*! This version of the function does not make use of the locator array. */
/*************************************************************************/
da_iakv_t da_iapqGetTopHeap(da_iapq_t *queue)
{
    ptr_t i, j;
    da_iakv_t *heap;
    idx_t node;
    val_t val;
    da_iakv_t vtx;

    if (queue->nnodes == 0)
        return (da_iakv_t){-1,0.0};

    queue->nnodes--;
    heap    = queue->heap;
    vtx.key = heap[0].key;
    vtx.val = heap[0].key;

    if ((i = queue->nnodes) > 0) {
        val  = heap[i].val;
        node = heap[i].key;
        i = 0;
        while ((j=2*i+1) < queue->nnodes) {
            if (val_lt(heap[j].val, val)) {
                if (j+1 < queue->nnodes && val_lt(heap[j+1].val, heap[j].val))
                    j = j+1;
                heap[i] = heap[j];
                i = j;
            }
            else if (j+1 < queue->nnodes && val_lt(heap[j+1].val, val)) {
                j = j+1;
                heap[i] = heap[j];
                i = j;
            }
            else
                break;
        }

        heap[i].val   = val;
        heap[i].key   = node;
    }

    return vtx;
}


/*************************************************************************/
/*! This function returns the item at the top of the queue. The item is not
    deleted from the queue. */
/**************************************************************************/
idx_t da_iapqSeeTopKey(da_iapq_t *queue)
{
    return (queue->nnodes == 0 ? -1 : queue->heap[0].key);
}


/*************************************************************************/
/*! This function returns the value of the top item. The item is not
    deleted from the queue. */
/**************************************************************************/
accum_t da_iapqSeeTopVal(da_iapq_t *queue)
{
    return (queue->nnodes == 0 ? 0.0 : queue->heap[0].val);
}


/*************************************************************************/
/*! This function returns the value of a specific item */
/**************************************************************************/
accum_t da_iapqSeeVal(da_iapq_t *queue, idx_t node)
{
    ptr_t *locator;
    da_iakv_t *heap;

    heap    = queue->heap;
    locator = queue->locator;

    return heap[locator[node]].val;
}

/*************************************************************************/
/*! This function returns the value of a specific item */
/**************************************************************************/
accum_t da_iapqSeeValHeap(da_iapq_t *queue, idx_t node)
{
    ptr_t i, j, loc;
    da_iakv_t *heap=queue->heap;

    /* linear search to find node */
    for(loc = -1, i=0; i < queue->nnodes; ++i)
        if(heap[i].key == node){
            loc = i;
            break;
        }

    return loc > -1 ? heap[loc].val : FLT_MAX;
}


/*************************************************************************/
/*! This function returns the value of a specific item */
/**************************************************************************/
int da_iapqExists(da_iapq_t *queue, idx_t node)
{
    ptr_t *locator;

    locator = queue->locator;

    return locator[node] != -1;
}

/*************************************************************************/
/*! This function returns the value of a specific item */
/**************************************************************************/
int da_iapqExistsHeap(da_iapq_t *queue, idx_t node)
{
    ptr_t i, j, loc;
    da_iakv_t *heap=queue->heap;

    /* linear search to find node */
    for(loc = -1, i=0; i < queue->nnodes; ++i)
        if(heap[i].key == node){
            loc = i;
            break;
        }

    return loc != -1;
}


/*************************************************************************/
/*! This functions checks the consistency of the heap */
/**************************************************************************/
int da_iapqCheckHeap(da_iapq_t *queue)
{
    ptr_t i, j;
    size_t nnodes;
    ptr_t *locator;
    da_iakv_t *heap;

    heap    = queue->heap;
    locator = queue->locator;
    nnodes  = queue->nnodes;

    if (!heap)
        da_errexit("CheckHeap: Heap not initialized!!!");

    if (nnodes == 0)
        return 1;

    ASSERT(locator[heap[0].key] == 0);
    for (i=1; i<nnodes; ++i) {
        ASSERT(locator[heap[i].key] == i);
        ASSERT(!val_lt(heap[i].val, heap[(i-1)/2].val));
        ASSERT(!val_lt(heap[i].val, heap[0].val));
    }

    for (j=i=0; i<queue->maxnodes; ++i) {
        if (locator[i] != -1)
            j++;
    }
    ASSERTP(j == nnodes, ("%jd %jd\n", (intmax_t)j, (intmax_t)nnodes));

    return 1;
}

void da_iapqPrintAll(da_iapq_t *queue){
    size_t i;

    printf("(%zu): ", queue->nnodes);
    for(i=0; i < queue->nnodes; ++i)
        printf(" %d %0.3f", queue->heap[i].key+1, queue->heap[i].val);
    printf("\n");

}


#undef val_lt
