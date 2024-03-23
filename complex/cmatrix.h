#ifndef _CMATRIX_
#define _CMATRIX_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

#define MIJ(SIZE, M, I, J) ((M)[(I) * (SIZE) + (J)])

#define copyAB(SIZE, A, B) do { \
	(void)memcpy((void *)(A), (void *)(B), (SIZE) * (SIZE) * sizeof (complex_t)); \
} while (0)

void mulCAB(int, complex_t *, complex_t *, complex_t *);
void mulCABtrans(int, complex_t *, complex_t *, complex_t *);
void mulCAtransB(int, complex_t *, complex_t *, complex_t *);
void transA(int, complex_t *);
void transAB(int, complex_t *, complex_t *);
void contransA(int, complex_t *); 
void contransAB(int, complex_t *, complex_t *);
void mulyAx(int, complex_t *, complex_t *, complex_t *);

#endif
