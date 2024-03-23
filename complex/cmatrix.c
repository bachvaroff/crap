#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cmatrix.h"

void
mulCAB(size, C, A, B)
	int size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	int i, j, k;

	for (i = 0u; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &MIJ(size, C, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(MIJ(size, C, i, j)) = 0.0;
			Im(MIJ(size, C, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(MIJ(size, C, i, j), MIJ(size, A, i, k), MIJ(size, B, k, j));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &MIJ(size, C, i, j));
			printcerr("", MIJ(size, C, i, j), "\n");
#endif
		}
	}
}

void
mulCAtransB(size, C, A, B)
	int size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	int i, j, k;

	for (i = 0u; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &MIJ(size, C, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(MIJ(size, C, i, j)) = 0.0;
			Im(MIJ(size, C, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(MIJ(size, C, i, j), MIJ(size, A, k, i), MIJ(size, B, k, j));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &MIJ(size, C, i, j));
			printcerr("", MIJ(size, C, i, j), "\n");
#endif
		}
	}
}

void
mulCABtrans(size, C, A, B)
	int size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	int i, j, k;

	for (i = 0u; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &MIJ(size, C, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(MIJ(size, C, i, j)) = 0.0;
			Im(MIJ(size, C, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(MIJ(size, C, i, j), MIJ(size, A, i, k), MIJ(size, B, j, k));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &MIJ(size, C, i, j));
			printcerr("", MIJ(size, C, i, j), "\n");
#endif
		}
	}
}

void
transA(size, A)
	int size;
	complex_t *A;
{
	int i, j;
	complex_t t;
	
	for (i = 0u; i < size; i++)
		for (j = i; j < size; j++) {
			t = MIJ(size, A, i, j);
			MIJ(size, A, i, j) = MIJ(size, A, j, i);
			MIJ(size, A, j, i) = t;
		}
	
	return;
}

void
transAB(size, A, B)
	int size;
	complex_t *A;
	complex_t *B;
{
	int i, j;
	
	for (i = 0u; i < size; i++)
		for (j = 0u; j < size; j++)
			MIJ(size, A, i, j) = MIJ(size, B, j, i);
	
	return;
}

void
contrans(size, A)
	int size;
	complex_t *A;
{
	int i, j;
	complex_t t;
	
	for (i = 0u; i < size; i++)
		for (j = i; j < size; j++) {
			con2(t, MIJ(size, A, i, j));
			con2(MIJ(size, A, i, j), MIJ(size, A, j, i));
			MIJ(size, A, j, i) = t;
		}
	
	return;
}

void
contransAB(size, A, B)
	int size;
	complex_t *A;
	complex_t *B;
{
	int i, j;
	
	for (i = 0u; i < size; i++)
		for (j = 0u; j < size; j++)
			con2(MIJ(size, A, i, j), MIJ(size, B, j, i));
	
	return;
}

void
mulYAX(size, Y, A, X)
	int size;
	complex_t *Y;
	complex_t *A;
	complex_t *X;
{
	int i, j;
	
	for (i = 0u; i < size; i++) {
		mkC(Y[i], 0.0, 0.0);
		for (j = 0u; j < size; j++)
			madd(Y[i], MIJ(size, A, i, j), X[j]);
	}
	
	return;
}
