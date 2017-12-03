/*!
\file  
\brief Various routines for providing portable 32 and 64 bit random number
       generators.

Ported from George Karypis' GKlib library by David C. Anastasiu, with permission, in Aug 2013.

\date   Started 5/17/2007
\author George
\version\verbatim $Id: random.c 14330 2013-05-18 12:15:15Z karypis $ \endverbatim
 */

#include "includes.h"

/*************************************************************************/
/*! GKlib's built in random number generator for portability across 
    different architectures */
/*************************************************************************/
#ifdef USE_DARAND
/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL	/* Most significant 33 bits */
#define LM 0x7FFFFFFFULL	/* Least significant 31 bits */


/* The array for the state vector */
static uint64_t mt[NN];
/* mti==NN+1 means mt[NN] is not initialized */
static int mti = NN + 1;
#endif /* USE_DARAND */

/* initializes mt[NN] with a seed */
void
da_randinit (uint64_t seed)
{
#ifdef USE_DARAND
  mt[0] = seed;
  for (mti = 1; mti < NN; mti++)
    mt[mti] =
      (6364136223846793005ULL * (mt[mti - 1] ^ (mt[mti - 1] >> 62)) + mti);
#endif
  srand ((unsigned int) seed);
}


/* generates a random number on [0, 2^64-1]-interval */
uint64_t
da_randint64 (void)
{
#ifdef USE_DARAND
  int i;
  unsigned long long x;
  static uint64_t mag01[2] = { 0ULL, MATRIX_A };

  if (mti >= NN)
    {				/* generate NN words at one time */
      /* if init_genrand64() has not been called, */
      /* a default initial seed is used     */
      if (mti == NN + 1)
	da_randinit (5489ULL);

      for (i = 0; i < NN - MM; i++)
	{
	  x = (mt[i] & UM) | (mt[i + 1] & LM);
	  mt[i] = mt[i + MM] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
	}
      for (; i < NN - 1; i++)
	{
	  x = (mt[i] & UM) | (mt[i + 1] & LM);
	  mt[i] = mt[i + (MM - NN)] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
	}
      x = (mt[NN - 1] & UM) | (mt[0] & LM);
      mt[NN - 1] = mt[MM - 1] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];

      mti = 0;
    }

  x = mt[mti++];

  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);

  return x & 0x7FFFFFFFFFFFFFFF;
#else
  return (uint64_t) (((uint64_t) rand ()) << 32 | ((uint64_t) rand ()));
#endif
}

/* generates a random number on [0, 2^32-1]-interval */
uint32_t
da_randint32 (void)
{
#ifdef USE_DARAND
  return (uint32_t) (da_randint64 () & 0x7FFFFFFF);
#else
  return (uint32_t) rand ();
#endif
}