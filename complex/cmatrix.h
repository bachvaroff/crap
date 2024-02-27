#ifndef _CMATRIX_
#define _CMATRIX_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

#define IDX(ARR, R, C, NC) ((ARR)[(R) * (NC) + (C)])

#define copyAB(SIZE, A, B) do { \
	(void)memcpy((void *)(A), (void *)(B), (SIZE) * (SIZE) * sizeof (complex_t)); \
} while (0)

void mulCAB(size_t, complex_t *, complex_t *, complex_t *);
void mulCABtrans(size_t, complex_t *, complex_t *, complex_t *);
void mulCAtransB(size_t, complex_t *, complex_t *, complex_t *);
void transA(size_t, complex_t *);
void transAB(size_t, complex_t *, complex_t *);
void contransA(size_t, complex_t *); 
void contransAB(size_t, complex_t *, complex_t *);

#endif
