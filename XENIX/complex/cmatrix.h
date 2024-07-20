#ifndef _CMATRIX_H
#define _CMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cconfig.h"
#include "complex.h"

#define MIJ(SIZE, M, I, J) ((M)[(I) * (SIZE) + (J)])

void printm(long, complex_t *, char *, int, int);
void printv(long, complex_t *, char *, int);
void printcv(long, complex_t *, char *, int);

/*
	cp - copy
	is - in situ
*/
								/* (N, A, B) : A <- B */
#define cpAB(SIZE, A, B) do { \
	(void)memcpy((void *)(A), (void *)(B), (SIZE) * (SIZE) * sizeof (complex_t)); \
} while (0)

void isAcA(long, complex_t *);					/* (N, A) : A <- conj(A) */
void cpAcB(long, complex_t *, complex_t *);			/* (N, A, B) : A <- conj(B) */ 
void isAtA(long, complex_t *);					/* (N, A) : A <- trans(A) */
void cpAtB(long, complex_t *, complex_t *);			/* (N, A, B) : A <- trans(B) */
void isActA(long, complex_t *);					/* (N, A) : A <- conj(trans(A)) */
void cpActB(long, complex_t *, complex_t *);			/* (N, A, B) : A <- conj(trans(B)) */
void mulCAB(long, complex_t *, complex_t *, complex_t *);	/* (N, C, A, B) : C <- A * B */
void mulCtAB(long, complex_t *, complex_t *, complex_t *);	/* (N, C, A, B) : C <- trans(A) * B */
void mulCctAB(long, complex_t *, complex_t *, complex_t *);	/* (N, C, A, B) : C <- conj(trans(A)) * B */
void mulCAtB(long, complex_t *, complex_t *, complex_t *);	/* (N, C, A, B) : C <- A * trans(B) */
void mulCActB(long, complex_t *, complex_t *, complex_t *);	/* (N, C, A, B) : C <- A * conj(trans(B)) */
void mulYvAXv(long, complex_t *, complex_t *, complex_t *);	/* (N, Yvec, A, Xvec) : Yvec <- A * Xvec */
void mulYcvXcvA(long, complex_t *, complex_t *, complex_t *);	/* (N, Ycovec, Xcovec, A) : Ycovec <- Xcovec * A */

#endif

