#ifndef _CMATRIX_H
#define _CMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

#define MIJ(SIZE, M, I, J) ((M)[(I) * (SIZE) + (J)])

#define copyAB(SIZE, A, B) do { \
	(void)memcpy((void *)(A), (void *)(B), (SIZE) * (SIZE) * sizeof (complex_t)); \
} while (0)

void mulCAB(long, complex_t *, complex_t *, complex_t *);
void mulCABtrans(long, complex_t *, complex_t *, complex_t *);
void mulCAtransB(long, complex_t *, complex_t *, complex_t *);
void mulCAcontransB(long, complex_t *, complex_t *, complex_t *);
void mulCABcontrans(long, complex_t *, complex_t *, complex_t *);
void transA(long, complex_t *);
void transAB(long, complex_t *, complex_t *);
void contransA(long, complex_t *); 
void contransAB(long, complex_t *, complex_t *);
void mulyAx(long, complex_t *, complex_t *, complex_t *);

#endif
