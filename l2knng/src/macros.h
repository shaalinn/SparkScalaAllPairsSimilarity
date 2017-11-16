/*
 * \file macros.h
 *
 *  Created on: Feb 3, 2014
 *      Author: dragos
 */

#ifndef _LIBL2AP_MACROS_H_
#define _LIBL2AP_MACROS_H_


/** MACROS  **/

#define F_NEQUALS_P(a,b,p) \
	do { if(fabs((a)-(b) ) > p) return 0; } while(0)


#define FL2INT(a,ll,f) \
	( (uint16_t) round((a - ll) * f) )

#define triget(m,i,j) \
	((i) >= (j) ? (m[((size_t)i)][((size_t)j)]) : (m[((size_t)j)][((size_t)i)]))

#define triset(m,i,j,v) \
	do { if((i)>=(j)){ m[((size_t)i)][((size_t)j)]=v; }else{ m[((size_t)j)][((size_t)i)]=v; } } while(0)

#define triincr(m,i,j) \
	do { if((i)>=(j)){ m[((size_t)i)][((size_t)j)]+=1; }else{ m[((size_t)j)][((size_t)i)]+=1; } } while(0)


#define COMPUTE_SCORE(f, idf, norm)\
    ( ((float)f)*idf/norm )



/********************
 * VMAT MACROS
 ********************/

/**
 * Notes: VMAT is a one-dimensional array, a condensed matrix of pairwise similarities or dissimilarities in
 * the same format as scipy.spatial.distance.pdist. It contains the flattened, upper-triangular part of a
 * pairwise (dis)similarity matrix, without the diagonal elements. If there are N data points and the matrix d
 * contains the (dis)similarity between the i-th and j-th observation at position d(i,j), the vector VMAT has
 * length N(N-1)/2 and is ordered as follows:
 *
 *   [ d(0,1), d(0,2), ..., d(0,n-1), d(1,2), ..., d(1,n-1), ..., d(n-2,n-1) ]
 */

/* get vmat[i][j] for vmat (upper triangular portion of dense matrix - diag) */
#define l2ap_vget(nrows, i, j) \
    ((i) >= (j) ? ( j*(2*nrows-j-3)/2+i-1 ) : ( i*(2*nrows-i-3)/2+j-1 ))
/* get nrows given length of vmat (condensed similarity matrix) */
#define l2ap_num_obs_y(vlen) \
    (int(ceil(sqrt(vlen * 2))))
/* get length of vmat (condensed similarity matrix) given nrows */
#define l2ap_vmat_len(nrows) \
    (nrows*(nrows-1)/2)


/*-------------------------------------------------------------
 * Usefull commands
 *-------------------------------------------------------------*/
#define da_max(a, b) ((a) >= (b) ? (a) : (b))
#define da_min(a, b) ((a) >= (b) ? (b) : (a))
#define da_max3(a, b, c) ((a) >= (b) && (a) >= (c) ? (a) : ((b) >= (a) && (b) >= (c) ? (b) : (c)))
#define DA_SWAP(a, b, tmp) do {(tmp) = (a); (a) = (b); (b) = (tmp);} while(0)
#define DA_INC_DEC(a, b, val) do {(a) += (val); (b) -= (val);} while(0)
#define da_sign(a, b) ((a >= 0 ? b : -b))

#define DA_ONEOVERRANDMAX (1.0/(RAND_MAX+1.0))
#define da_randomInRange(u) ((int) (DA_ONEOVERRANDMAX*(u)*rand()))
#define da_andomInRange_r(s, u) ((int) (DA_ONEOVERRANDMAX*(u)*rand_r(s)))

#define da_abs(x) ((x) >= 0 ? (x) : -(x))


/*-------------------------------------------------------------
 * Timing macros
 *-------------------------------------------------------------*/
#define da_clearcputimer(tmr) (tmr = 0.0)
#define da_startcputimer(tmr) (tmr -= da_CPUSeconds())
#define da_stopcputimer(tmr)  (tmr += da_CPUSeconds())
#define da_getcputimer(tmr)   (tmr)

#define da_clearwctimer(tmr) (tmr = 0.0)
#define da_startwctimer(tmr) (tmr -= da_WClockSeconds())
#define da_stopwctimer(tmr)  (tmr += da_WClockSeconds())
#define da_getwctimer(tmr)   (tmr)


/**
 * OMP macros
 */

/**
 * \param myid thread id
 * \param val private var doing max over
 * \param buffer of same type as val and size numThreads
 * \param nthreads number of threads
 */
#define da_omp_maxreduce(myid,val,buffer,nthreads) \
  do { \
    _Pragma("omp barrier") \
    int32_t _i; \
    (buffer)[myid] = (val); \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      for (_i=0;_i<(nthreads);++_i) { \
        if ((val) < (buffer)[_i]) { \
          (val) = (buffer)[_i]; \
        } \
      } \
      (buffer)[0] = (val); \
    } \
    _Pragma("omp barrier") \
    (val) = (buffer)[0]; \
  } while(0)

#define da_omp_minreduce(myid,val,buffer,nthreads) \
  do { \
    _Pragma("omp barrier") \
    int32_t _i; \
    (buffer)[myid] = (val); \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      for (_i=0;_i<(nthreads);++_i) { \
        if ((val) < (buffer)[_i]) { \
          (val) = (buffer)[_i]; \
        } \
      } \
      (buffer)[0] = (val); \
    } \
    _Pragma("omp barrier") \
    (val) = (buffer)[0]; \
  } while(0)



/*-------------------------------------------------------------
 * CSR conversion macros
 *-------------------------------------------------------------*/
#define CSRMAKE(i, n, a) \
   do { \
     for (i=1; i<n; i++) a[i] += a[i-1]; \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0)

#define CSRSHIFT(i, n, a) \
   do { \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0)

/*********************
 * Progress indicator
 *********************/

#define NPCT    10   /* number of steps */
#define da_progress_init(pct, indicator, niterations) \
    do { \
        indicator = ceil(niterations/(float)NPCT); \
        pct = 0; \
    } while(0)

#define da_progress_init_steps(pct, indicator, niterations, nsteps) \
    do { \
        indicator = ceil(niterations/(float)nsteps); \
        pct = 0; \
    } while(0)

#define da_progress_advance(pct) \
    do { \
        if(pct > 0 && pct < 100) \
            printf("%d%%..", pct); \
        fflush(stdout); \
        pct += NPCT; \
    } while(0)

#define da_progress_advance_steps(pct, nsteps) \
    do { \
        if(pct > 0 && pct < 100) \
            printf("%d%%..", pct); \
        fflush(stdout); \
        pct += 100.0/nsteps; \
    } while(0)

#define da_progress_finalize(pct) \
    do { \
        while(pct < 100){ \
            if(pct > 0) \
                printf("%d%%..", pct); \
            pct += NPCT; \
        } \
        if(pct == 100) \
            printf("%d%%", pct); \
        fflush(stdout); \
    } while(0)

#define da_progress_finalize_steps(pct, nsteps) \
    do { \
        while(pct < 100){ \
            if(pct > 0) \
                printf("%d%%..", pct); \
            pct += 100.0/nsteps; \
        } \
        if(pct == 100) \
            printf("%d%%", pct); \
        fflush(stdout); \
    } while(0)

/**
 * Debug macros
 */

#ifndef NDEBUG
#   define ASSERT(expr)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        assert(expr);                                                \
    }

#   define ASSERTP(expr,msg)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        printf msg ; \
        printf("\n"); \
        assert(expr);                                                \
    }
#else
#   define ASSERT(expr) ;
#   define ASSERTP(expr,msg) ;
#endif


#define __BACKTRACE() \
  do { \
    void * _buffer[255]; \
    int _size = backtrace(_buffer,255); \
    char ** _strings = backtrace_symbols(_buffer,_size); \
    int _i; \
    fprintf(stderr,"\n"); \
    for (_i=0;_i<_size;++_i) { \
      fprintf(stderr,"%d:[%p] %s\n",_i,_buffer[_i],_strings[_i]); \
    } \
    free(_strings); \
  } while(0)

#define eprintf(fmt, ...) \
  do { \
    fprintf(stderr , "ERROR: "fmt , ##__VA_ARGS__); \
    fflush(stderr); \
  } while (0)

#ifndef NDEBUG
  #define wprintf(fmt, ...) \
    do { \
      fprintf(stderr , "WARN: "fmt , ##__VA_ARGS__); \
      fflush(stderr); \
    } while (0)

  #define dprintf(fmt, ...) \
    do { \
      fprintf(stdout , "DEBUG: "fmt , ##__VA_ARGS__); \
      fflush(stdout); \
    } while (0)

  #define DA_ASSERT(cond,fmt, ...) \
    do { \
      if (!(cond)) { \
        eprintf(fmt, ##__VA_ARGS__); \
        /* print the stack trace */ \
        __BACKTRACE(); \
        assert(0); /* will always fail */ \
      } \
    } while (0)

  #define DA_ASSERT_EQUALS(a,b,fmt) \
    DA_ASSERT(a == b,"("#a" = "fmt") != ("#b" = "fmt")",a,b)
#else
	#define wprintf(fmt, ...)
	#define dprintf(fmt, ...)
	#define DL_ASSERT(cnd,fmt, ...)
	#define DL_ASSERT_EQUALS(a,b,fmt)
#endif



#endif /* MACROS_H_ */
