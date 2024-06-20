#ifndef _CLU_H
#define _CLU_H

#include "cconfig.h"
#include "complex.h"
#include "cmatrix.h"

int LUPDecompose(long, complex_t *, long *, double);
void LUPExtract(long, complex_t *, complex_t *, complex_t *,
	complex_t *, long *);
void LUPSolve(long, complex_t *, long *, complex_t *, complex_t *);
void LUPInvert(long, complex_t *, complex_t *, long *);
complex_t LUPDeterminant(long, complex_t *, long *);

#endif

