#ifndef _CLU_
#define _CLU_

#include "cmatrix.h"

int LUPDecompose(long, complex_t *, long *, double);
void LUPSolve(long, complex_t *, long *, complex_t *, complex_t *);
void LUPInvert(long, complex_t *, complex_t *, long *);
complex_t LUPDeterminant(long, complex_t *, long *);

#endif

