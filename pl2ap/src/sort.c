/*!
\file  sort.c
\brief This file contains various sorting routines

These routines are implemented using the DASORT macro that is defined
in da_qsort.h and is based on GNU's GLIBC qsort() implementation.

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Aug 2013.

Additional sorting routines can be created using the same way that
these routines where defined.

\date   Started 2/27/2013
\author David C. Anastasiu
 */

#include "includes.h"



/*************************************************************************/
/*! Sorts an array of ptr_t in increasing order */
/*************************************************************************/
void
da_psorti (const size_t n, ptr_t * const base)
{
#define ptr_lt(a, b) ((*a) < (*b))
  DA_MKQSORT (ptr_t, base, n, ptr_lt);
#undef ptr_lt
}


/*************************************************************************/
/*! Sorts an array of ptr_t in decreasing order */
/*************************************************************************/
void
da_psortd (const size_t n, ptr_t * const base)
{
#define ptr_gt(a, b) ((*a) > (*b))
  DA_MKQSORT (ptr_t, base, n, ptr_gt);
#undef ptr_gt
}


/*************************************************************************/
/*! Sorts an array of idx_t in increasing order */
/*************************************************************************/
void
da_isorti (const size_t n, idx_t * const base)
{
#define idx_lt(a, b) ((*a) < (*b))
  DA_MKQSORT (idx_t, base, n, idx_lt);
#undef idx_lt
}


/*************************************************************************/
/*! Sorts an array of idx_t in decreasing order */
/*************************************************************************/
void
da_isortd (const size_t n, idx_t * const base)
{
#define idx_gt(a, b) ((*a) > (*b))
  DA_MKQSORT (idx_t, base, n, idx_gt);
#undef idx_gt
}


/*************************************************************************/
/*! Sorts an array of val_t in increasing order */
/*************************************************************************/
void
da_vsorti (const size_t n, val_t * const base)
{
#define val_lt(a, b) ((*a) < (*b))
  DA_MKQSORT (val_t, base, n, val_lt);
#undef val_lt
}


/*************************************************************************/
/*! Sorts an array of val_t in decreasing order */
/*************************************************************************/
void
da_vsortd (const size_t n, val_t * const base)
{
#define float_gt(a, b) ((*a) > (*b))
  DA_MKQSORT (val_t, base, n, float_gt);
#undef float_gt
}



/*************************************************************************/
/*! Sorts an array of da_ppkv_t in increasing order */
/*************************************************************************/
void
da_ppkvsorti (const size_t n, da_ppkv_t * const base)
{
#define da_ppkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ppkv_t, base, n, da_ppkv_lt);
#undef da_ppkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ppkv_t in decreasing order */
/*************************************************************************/
void
da_ppkvsortd (const size_t n, da_ppkv_t * const base)
{
#define da_ppkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ppkv_t, base, n, da_ppkv_gt);
#undef da_ppkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pikv_t in increasing order */
/*************************************************************************/
void
da_pikvsorti (const size_t n, da_pikv_t * const base)
{
#define da_pikv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pikv_t, base, n, da_pikv_lt);
#undef da_pikv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pikv_t in decreasing order */
/*************************************************************************/
void
da_pikvsortd (const size_t n, da_pikv_t * const base)
{
#define da_pikv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pikv_t, base, n, da_pikv_gt);
#undef da_pikv_gt
}



/*************************************************************************/
/*! Sorts an array of da_pvkv_t in increasing order */
/*************************************************************************/
void
da_pvkvsorti (const size_t n, da_pvkv_t * const base)
{
#define da_pvkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pvkv_t, base, n, da_pvkv_lt);
#undef da_pvkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pvkv_t in decreasing order */
/*************************************************************************/
void
da_pvkvsortd (const size_t n, da_pvkv_t * const base)
{
#define da_pvkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pvkv_t, base, n, da_pvkv_gt);
#undef da_pvkv_gt
}


/*************************************************************************/
/*! Sorts an array of da_pckv_t in increasing order */
/*************************************************************************/
void
da_pckvsorti (const size_t n, da_pckv_t * const base)
{
#define da_pckv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pckv_t, base, n, da_pckv_lt);
#undef da_pckv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pckv_t in decreasing order */
/*************************************************************************/
void
da_pckvsortd (const size_t n, da_pckv_t * const base)
{
#define da_pckv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pckv_t, base, n, da_pckv_gt);
#undef da_pkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pi32kv_t in increasing order */
/*************************************************************************/
void
da_pi32kvsorti (const size_t n, da_pi32kv_t * const base)
{
#define da_pi32kv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pi32kv_t, base, n, da_pi32kv_lt);
#undef da_pi32kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pi32kv_t in decreasing order */
/*************************************************************************/
void
da_pi32kvsortd (const size_t n, da_pi32kv_t * const base)
{
#define da_pi32kv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pi32kv_t, base, n, da_pi32kv_gt);
#undef da_pi32kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pi64kv_t in increasing order */
/*************************************************************************/
void
da_pi64kvsorti (const size_t n, da_pi64kv_t * const base)
{
#define da_pi64kv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pi64kv_t, base, n, da_pi64kv_lt);
#undef da_pi64kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pi64kv_t in decreasing order */
/*************************************************************************/
void
da_pi64kvsortd (const size_t n, da_pi64kv_t * const base)
{
#define da_pi64kv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pi64kv_t, base, n, da_pi64kv_gt);
#undef da_pi64kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pzkv_t in increasing order */
/*************************************************************************/
void
da_pzkvsorti (const size_t n, da_pzkv_t * const base)
{
#define da_pzkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pzkv_t, base, n, da_pzkv_lt);
#undef da_pzkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pzkv_t in decreasing order */
/*************************************************************************/
void
da_pzkvsortd (const size_t n, da_pzkv_t * const base)
{
#define da_pzkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pzkv_t, base, n, da_pzkv_gt);
#undef da_pzkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pfkv_t in increasing order */
/*************************************************************************/
void
da_pfkvsorti (const size_t n, da_pfkv_t * const base)
{
#define da_pfkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pfkv_t, base, n, da_pfkv_lt);
#undef da_pfkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pfkv_t in decreasing order */
/*************************************************************************/
void
da_pfkvsortd (const size_t n, da_pfkv_t * const base)
{
#define da_pfkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pfkv_t, base, n, da_pfkv_gt);
#undef da_pfkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pdkv_t in increasing order */
/*************************************************************************/
void
da_pdkvsorti (const size_t n, da_pdkv_t * const base)
{
#define da_pdkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pdkv_t, base, n, da_pdkv_lt);
#undef da_pdkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pdkv_t in decreasing order */
/*************************************************************************/
void
da_pdkvsortd (const size_t n, da_pdkv_t * const base)
{
#define da_pdkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pdkv_t, base, n, da_pdkv_gt);
#undef da_pdkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pskv_t in increasing order */
/*************************************************************************/
void
da_pskvsorti (const size_t n, da_pskv_t * const base)
{
#define da_pskv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_pskv_t, base, n, da_pskv_lt);
#undef da_pskv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pskv_t in decreasing order */
/*************************************************************************/
void
da_pskvsortd (const size_t n, da_pskv_t * const base)
{
#define da_pskv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_pskv_t, base, n, da_pskv_gt);
#undef da_pskv_gt
}


/* da_i ## kv pairs */


/*************************************************************************/
/*! Sorts an array of da_ipkv_t in increasing order */
/*************************************************************************/
void
da_ipkvsorti (const size_t n, da_ipkv_t * const base)
{
#define da_ipkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ipkv_t, base, n, da_ipkv_lt);
#undef da_ipkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ipkv_t in decreasing order */
/*************************************************************************/
void
da_ipkvsortd (const size_t n, da_ipkv_t * const base)
{
#define da_ipkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ipkv_t, base, n, da_ipkv_gt);
#undef da_ipkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pikv_t in increasing order */
/*************************************************************************/
void
da_iikvsorti (const size_t n, da_iikv_t * const base)
{
#define da_iikv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_iikv_t, base, n, da_iikv_lt);
#undef da_iikv_lt
}


/*************************************************************************/
/*! Sorts an array of da_iikv_t in decreasing order */
/*************************************************************************/
void
da_iikvsortd (const size_t n, da_iikv_t * const base)
{
#define da_iikv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_iikv_t, base, n, da_iikv_gt);
#undef da_iikv_gt
}



/*************************************************************************/
/*! Sorts an array of da_ivkv_t in increasing order */
/*************************************************************************/
void
da_ivkvsorti (const size_t n, da_ivkv_t * const base)
{
#define da_ivkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ivkv_t, base, n, da_ivkv_lt);
#undef da_ivkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ivkv_t in decreasing order */
/*************************************************************************/
void
da_ivkvsortd (const size_t n, da_ivkv_t * const base)
{
#define da_ivkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ivkv_t, base, n, da_ivkv_gt);
#undef da_ivkv_gt
}



/*************************************************************************/
/*! Sorts an array of da_iakv_t in increasing order */
/*************************************************************************/
void
da_iakvsorti (const size_t n, da_iakv_t * const base)
{
#define da_iakv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_iakv_t, base, n, da_iakv_lt);
#undef da_iakv_lt
}


/*************************************************************************/
/*! Sorts an array of da_iakv_t in decreasing order */
/*************************************************************************/
void
da_iakvsortd (const size_t n, da_iakv_t * const base)
{
#define da_iakv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_iakv_t, base, n, da_iakv_gt);
#undef da_iakv_gt
}





/*************************************************************************/
/*! Sorts an array of da_ickv_t in increasing order */
/*************************************************************************/
void
da_ickvsorti (const size_t n, da_ickv_t * const base)
{
#define da_ickv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ickv_t, base, n, da_ickv_lt);
#undef da_ickv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ickv_t in decreasing order */
/*************************************************************************/
void
da_ickvsortd (const size_t n, da_ickv_t * const base)
{
#define da_ickv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ickv_t, base, n, da_ickv_gt);
#undef da_pkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ii32kv_t in increasing order */
/*************************************************************************/
void
da_ii32kvsorti (const size_t n, da_ii32kv_t * const base)
{
#define da_ii32kv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ii32kv_t, base, n, da_ii32kv_lt);
#undef da_ii32kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ii32kv_t in decreasing order */
/*************************************************************************/
void
da_ii32kvsortd (const size_t n, da_ii32kv_t * const base)
{
#define da_ii32kv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ii32kv_t, base, n, da_ii32kv_gt);
#undef da_ii32kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ii64kv_t in increasing order */
/*************************************************************************/
void
da_ii64kvsorti (const size_t n, da_ii64kv_t * const base)
{
#define da_ii64kv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ii64kv_t, base, n, da_ii64kv_lt);
#undef da_ii64kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ii64kv_t in decreasing order */
/*************************************************************************/
void
da_ii64kvsortd (const size_t n, da_ii64kv_t * const base)
{
#define da_ii64kv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ii64kv_t, base, n, da_ii64kv_gt);
#undef da_ii64kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_izkv_t in increasing order */
/*************************************************************************/
void
da_izkvsorti (const size_t n, da_izkv_t * const base)
{
#define da_izkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_izkv_t, base, n, da_izkv_lt);
#undef da_izkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_izkv_t in decreasing order */
/*************************************************************************/
void
da_izkvsortd (const size_t n, da_izkv_t * const base)
{
#define da_izkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_izkv_t, base, n, da_izkv_gt);
#undef da_izkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ifkv_t in increasing order */
/*************************************************************************/
void
da_ifkvsorti (const size_t n, da_ifkv_t * const base)
{
#define da_ifkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ifkv_t, base, n, da_ifkv_lt);
#undef da_ifkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ifkv_t in decreasing order */
/*************************************************************************/
void
da_ifkvsortd (const size_t n, da_ifkv_t * const base)
{
#define da_ifkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ifkv_t, base, n, da_ifkv_gt);
#undef da_ifkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_idkv_t in increasing order */
/*************************************************************************/
void
da_idkvsorti (const size_t n, da_idkv_t * const base)
{
#define da_idkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_idkv_t, base, n, da_idkv_lt);
#undef da_idkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_idkv_t in decreasing order */
/*************************************************************************/
void
da_idkvsortd (const size_t n, da_idkv_t * const base)
{
#define da_idkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_idkv_t, base, n, da_idkv_gt);
#undef da_idkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_iskv_t in increasing order */
/*************************************************************************/
void
da_iskvsorti (const size_t n, da_iskv_t * const base)
{
#define da_iskv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_iskv_t, base, n, da_iskv_lt);
#undef da_iskv_lt
}


/*************************************************************************/
/*! Sorts an array of da_iskv_t in decreasing order */
/*************************************************************************/
void
da_iskvsortd (const size_t n, da_iskv_t * const base)
{
#define da_iskv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_iskv_t, base, n, da_iskv_gt);
#undef da_iskv_gt
}






/* da_u ## kv pairs */


/*************************************************************************/
/*! Sorts an array of da_upkv_t in increasing order */
/*************************************************************************/
void
da_upkvsorti (const size_t n, da_upkv_t * const base)
{
#define da_upkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_upkv_t, base, n, da_upkv_lt);
#undef da_upkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_upkv_t in decreasing order */
/*************************************************************************/
void
da_upkvsortd (const size_t n, da_upkv_t * const base)
{
#define da_upkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_upkv_t, base, n, da_upkv_gt);
#undef da_upkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pikv_t in increasing order */
/*************************************************************************/
void
da_uikvsorti (const size_t n, da_uikv_t * const base)
{
#define da_uikv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_uikv_t, base, n, da_uikv_lt);
#undef da_uikv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uikv_t in decreasing order */
/*************************************************************************/
void
da_uikvsortd (const size_t n, da_uikv_t * const base)
{
#define da_uikv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_uikv_t, base, n, da_uikv_gt);
#undef da_uikv_gt
}



/*************************************************************************/
/*! Sorts an array of da_uvkv_t in increasing order */
/*************************************************************************/
void
da_uvkvsorti (const size_t n, da_uvkv_t * const base)
{
#define da_uvkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_uvkv_t, base, n, da_uvkv_lt);
#undef da_uvkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uvkv_t in decreasing order */
/*************************************************************************/
void
da_uvkvsortd (const size_t n, da_uvkv_t * const base)
{
#define da_uvkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_uvkv_t, base, n, da_uvkv_gt);
#undef da_uvkv_gt
}



/*************************************************************************/
/*! Sorts an array of da_uakv_t in increasing order */
/*************************************************************************/
void
da_uakvsorti (const size_t n, da_uakv_t * const base)
{
#define da_uakv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_uakv_t, base, n, da_uakv_lt);
#undef da_uakv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uakv_t in decreasing order */
/*************************************************************************/
void
da_uakvsortd (const size_t n, da_uakv_t * const base)
{
#define da_uakv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_uakv_t, base, n, da_uakv_gt);
#undef da_uakv_gt
}



/*************************************************************************/
/*! Sorts an array of da_uckv_t in increasing order */
/*************************************************************************/
void
da_uckvsorti (const size_t n, da_uckv_t * const base)
{
#define da_uckv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_uckv_t, base, n, da_uckv_lt);
#undef da_uckv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uckv_t in decreasing order */
/*************************************************************************/
void
da_uckvsortd (const size_t n, da_uckv_t * const base)
{
#define da_uckv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_uckv_t, base, n, da_uckv_gt);
#undef da_pkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ui32kv_t in increasing order */
/*************************************************************************/
void
da_ui32kvsorti (const size_t n, da_ui32kv_t * const base)
{
#define da_ui32kv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ui32kv_t, base, n, da_ui32kv_lt);
#undef da_ui32kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ui32kv_t in decreasing order */
/*************************************************************************/
void
da_ui32kvsortd (const size_t n, da_ui32kv_t * const base)
{
#define da_ui32kv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ui32kv_t, base, n, da_ui32kv_gt);
#undef da_ui32kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ui64kv_t in increasing order */
/*************************************************************************/
void
da_ui64kvsorti (const size_t n, da_ui64kv_t * const base)
{
#define da_ui64kv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ui64kv_t, base, n, da_ui64kv_lt);
#undef da_ui64kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ui64kv_t in decreasing order */
/*************************************************************************/
void
da_ui64kvsortd (const size_t n, da_ui64kv_t * const base)
{
#define da_ui64kv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ui64kv_t, base, n, da_ui64kv_gt);
#undef da_ui64kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_uzkv_t in increasing order */
/*************************************************************************/
void
da_uzkvsorti (const size_t n, da_uzkv_t * const base)
{
#define da_uzkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_uzkv_t, base, n, da_uzkv_lt);
#undef da_uzkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uzkv_t in decreasing order */
/*************************************************************************/
void
da_uzkvsortd (const size_t n, da_uzkv_t * const base)
{
#define da_uzkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_uzkv_t, base, n, da_uzkv_gt);
#undef da_uzkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ufkv_t in increasing order */
/*************************************************************************/
void
da_ufkvsorti (const size_t n, da_ufkv_t * const base)
{
#define da_ufkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_ufkv_t, base, n, da_ufkv_lt);
#undef da_ufkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ufkv_t in decreasing order */
/*************************************************************************/
void
da_ufkvsortd (const size_t n, da_ufkv_t * const base)
{
#define da_ufkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_ufkv_t, base, n, da_ufkv_gt);
#undef da_ufkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_udkv_t in increasing order */
/*************************************************************************/
void
da_udkvsorti (const size_t n, da_udkv_t * const base)
{
#define da_udkv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_udkv_t, base, n, da_udkv_lt);
#undef da_udkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_udkv_t in decreasing order */
/*************************************************************************/
void
da_udkvsortd (const size_t n, da_udkv_t * const base)
{
#define da_udkv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_udkv_t, base, n, da_udkv_gt);
#undef da_udkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_uskv_t in increasing order */
/*************************************************************************/
void
da_uskvsorti (const size_t n, da_uskv_t * const base)
{
#define da_uskv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_uskv_t, base, n, da_uskv_lt);
#undef da_uskv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uskv_t in decreasing order */
/*************************************************************************/
void
da_uskvsortd (const size_t n, da_uskv_t * const base)
{
#define da_uskv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_uskv_t, base, n, da_uskv_gt);
#undef da_uskv_gt
}



/*************************************************************************/
/*! Sorts an array of da_il_t in increasing order */
/*************************************************************************/
void
da_ilsorti (const size_t n, da_il_t ** const base)
{
#define da_il_lt(a, b) ((*a)->maxscore < (*b)->maxscore)
  DA_MKQSORT (da_il_t *, base, n, da_il_lt);
#undef da_il_lt
}


/*************************************************************************/
/*! Sorts an array of da_il_t in decreasing order */
/*************************************************************************/
void
da_ilsortd (const size_t n, da_il_t ** const base)
{
#define da_il_gt(a, b) ((*a)->maxscore > (*b)->maxscore)
  DA_MKQSORT (da_il_t *, base, n, da_il_gt);
#undef da_il_gt
}



/*************************************************************************/
/*! Sorts an array of da_bil_t in increasing order */
/*************************************************************************/
void
da_bilsorti (const size_t n, da_bil_t ** const base)
{
#define da_bil_lt(a, b) ((*a)->maxscore < (*b)->maxscore)
  DA_MKQSORT (da_bil_t *, base, n, da_bil_lt);
#undef da_bil_lt
}


/*************************************************************************/
/*! Sorts an array of da_bil_t in decreasing order */
/*************************************************************************/
void
da_bilsortd (const size_t n, da_bil_t ** const base)
{
#define da_bil_gt(a, b) ((*a)->maxscore > (*b)->maxscore)
  DA_MKQSORT (da_bil_t *, base, n, da_bil_gt);
#undef da_bil_gt
}



/*************************************************************************/
/*! Sorts an array of da_uukv_t in increasing order */
/*************************************************************************/
void
da_uukvsorti (const size_t n, da_uukv_t * const base)
{
#define da_uukv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_uukv_t, base, n, da_uukv_lt);
#undef da_uukv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uukv_t in decreasing order */
/*************************************************************************/
void
da_uukvsortd (const size_t n, da_uukv_t * const base)
{
#define da_uukv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_uukv_t, base, n, da_uukv_gt);
#undef da_uukv_gt
}




/*************************************************************************/
/*! Sorts an array of da_iiakv_t in increasing order */
/*************************************************************************/
void
da_iiakvsorti (const size_t n, da_iiakv_t * const base)
{
#define da_iiakv_lt(a, b) ((a)->val < (b)->val)
  DA_MKQSORT (da_iiakv_t, base, n, da_iiakv_lt);
#undef da_iiakv_lt
}


/*************************************************************************/
/*! Sorts an array of da_iiakv_t in decreasing order */
/*************************************************************************/
void
da_iiakvsortd (const size_t n, da_iiakv_t * const base)
{
#define da_iiakv_gt(a, b) ((a)->val > (b)->val)
  DA_MKQSORT (da_iiakv_t, base, n, da_iiakv_gt);
#undef da_iiakv_gt
}


/*************************************************************************/
/*! Sorts an array of da_iiakv_t in increasing order by keys */
/*************************************************************************/
void
da_iiakeysorti (const size_t n, da_iiakv_t * const base)
{
#define da_iiakv_lt(a, b) ((a)->key1 < (b)->key1 || ((a)->key1 == (b)->key1 && (a)->key2 < (b)->key2 ))
  DA_MKQSORT (da_iiakv_t, base, n, da_iiakv_lt);
#undef da_iiakv_lt
}


/*************************************************************************/
/*! Sorts an array of da_iiakv_t in decreasing order by keys */
/*************************************************************************/
void
da_iiakeysortd (const size_t n, da_iiakv_t * const base)
{
#define da_iiakv_gt(a, b) ((a)->key1 > (b)->key1 || ((a)->key1 == (b)->key1 && (a)->key2 > (b)->key2 ))
  DA_MKQSORT (da_iiakv_t, base, n, da_iiakv_gt);
#undef da_iiakv_gt
}
