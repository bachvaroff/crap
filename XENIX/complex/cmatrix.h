#ifndef _CMATRIX_H
#define _CMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

#define MIJ(SIZE, M, I, J) ((M)[(I) * (SIZE) + (J)])

									/* (N, A, B) : A <- B */
#define AcopyB(SIZE, A, B) do { \
	(void)memcpy((void *)(A), (void *)(B), (SIZE) * (SIZE) * sizeof (complex_t)); \
} while (0)

void AconjA(long, complex_t *);						/* (N, A) : A <- conj(A) */
void AconjB(long, complex_t *, complex_t *);				/* (N, A, B) : A <- conj(B) */ 
void AtransA(long, complex_t *);					/* (N, A) : A <- trans(A) */
void AtransB(long, complex_t *, complex_t *);				/* (N, A, B) : A <- trans(B) */
void AconjtransA(long, complex_t *);					/* (N, A) : A <- conj(trans(A)) */
void AconjtransB(long, complex_t *, complex_t *);			/* (N, A, B) : A <- conj(trans(B)) */
void mulCAB(long, complex_t *, complex_t *, complex_t *);		/* (N, C, A, B) : C <- A * B */
void mulCAtransB(long, complex_t *, complex_t *, complex_t *);		/* (N, C, A, B) : C <- A * trans(B) */
void mulCtransAB(long, complex_t *, complex_t *, complex_t *);		/* (N, C, A, B) : C <- trans(A) * B */
void mulCconjtransAB(long, complex_t *, complex_t *, complex_t *);	/* (N, C, A, B) : C <- conj(trans(A)) * B */
void mulCAconjtransB(long, complex_t *, complex_t *, complex_t *);	/* (N, C, A, B) : C <- A * conj(trans(B)) */
void mulyvAxv(long, complex_t *, complex_t *, complex_t *);		/* (N, yvec, A, xvec) : yvec <- A * xvec */
void mulycovxcovA(long, complex_t *, complex_t *, complex_t *);		/* (N, ycovec, xcovec, A) : ycovec <- xcovec * A */

#endif
