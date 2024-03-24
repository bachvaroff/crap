#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cmatrix.h"

void
mulCAB(size, C, A, B)
	long size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	long i, j, k;

	for (i = 0l; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%16lx % 11ld %16lx\n",
			(unsigned long)C, i,
			(unsigned long)&MIJ(size, C, i, 0l));
#endif
		for (j = 0l; j < size; j++) {
			Re(MIJ(size, C, i, j)) = 0.0;
			Im(MIJ(size, C, i, j)) = 0.0;
			for (k = 0l; k < size; k++)
				madd(MIJ(size, C, i, j), MIJ(size, A, i, k), MIJ(size, B, k, j));
		}
	}
}

void
mulCAtransB(size, C, A, B)
	long size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	long i, j, k;

	for (i = 0l; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%16lx % 11ld %16lx\n",
			(unsigned long)C, i,
			(unsigned long)&MIJ(size, C, i, 0l));
#endif
		for (j = 0l; j < size; j++) {
			Re(MIJ(size, C, i, j)) = 0.0;
			Im(MIJ(size, C, i, j)) = 0.0;
			for (k = 0l; k < size; k++)
				madd(MIJ(size, C, i, j), MIJ(size, A, k, i), MIJ(size, B, k, j));
		}
	}
}

void
mulCABtrans(size, C, A, B)
	long size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	long i, j, k;

	for (i = 0l; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%16lx % 11ld %16lx\n",
			(unsigned long)C, i,
			(unsigned long)&MIJ(size, C, i, 0l));
#endif
		for (j = 0l; j < size; j++) {
			Re(MIJ(size, C, i, j)) = 0.0;
			Im(MIJ(size, C, i, j)) = 0.0;
			for (k = 0l; k < size; k++)
				madd(MIJ(size, C, i, j), MIJ(size, A, i, k), MIJ(size, B, j, k));
		}
	}
}

void
transA(size, A)
	long size;
	complex_t *A;
{
	long i, j;
	complex_t t;
	
	for (i = 0l; i < size; i++)
		for (j = i; j < size; j++) {
			t = MIJ(size, A, i, j);
			MIJ(size, A, i, j) = MIJ(size, A, j, i);
			MIJ(size, A, j, i) = t;
		}
	
	return;
}

void
transAB(size, A, B)
	long size;
	complex_t *A;
	complex_t *B;
{
	long i, j;
	
	for (i = 0l; i < size; i++)
		for (j = 0l; j < size; j++)
			MIJ(size, A, i, j) = MIJ(size, B, j, i);
	
	return;
}

void
contrans(size, A)
	long size;
	complex_t *A;
{
	long i, j;
	complex_t t;
	
	for (i = 0l; i < size; i++)
		for (j = i; j < size; j++) {
			con2(t, MIJ(size, A, i, j));
			con2(MIJ(size, A, i, j), MIJ(size, A, j, i));
			MIJ(size, A, j, i) = t;
		}
	
	return;
}

void
contransAB(size, A, B)
	long size;
	complex_t *A;
	complex_t *B;
{
	long i, j;
	
	for (i = 0l; i < size; i++)
		for (j = 0l; j < size; j++)
			con2(MIJ(size, A, i, j), MIJ(size, B, j, i));
	
	return;
}

void
mulyAx(size, y, A, x)
	long size;
	complex_t *y;
	complex_t *A;
	complex_t *x;
{
	long i, j;
	
	for (i = 0l; i < size; i++) {
		mkC(y[i], 0.0, 0.0);
		for (j = 0l; j < size; j++)
			madd(y[i], MIJ(size, A, i, j), x[j]);
	}
	
	return;
}

