#ifndef _CLU_
#define _CLU_

#include "cmatrix.h"

int LUPDecompose(int, complex_t *, double, int *);
void LUPSolve(int, complex_t *, int *, complex_t *, complex_t *);
void LUPInvert(int, complex_t *, int *, complex_t *);
complex_t LUPDeterminant(int, complex_t *, int *);

#endif

