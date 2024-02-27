#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cmatrix.h"

void
mulCAB(size, C, A, B)
	size_t size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	size_t i, j, k;

	for (i = 0u; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &IDX(C, size, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(IDX(C, size, i, j)) = 0.0;
			Im(IDX(C, size, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(IDX(C, size, i, j), IDX(A, size, i, k), IDX(B, size, k, j));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &IDX(C, size, i, j));
			printcerr("", IDX(C, size, i, j), "\n");
#endif
		}
	}
}

void
mulCAtransB(size, C, A, B)
	size_t size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	size_t i, j, k;

	for (i = 0u; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &IDX(C, size, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(IDX(C, size, i, j)) = 0.0;
			Im(IDX(C, size, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(IDX(C, size, i, j), IDX(A, size, k, i), IDX(B, size, k, j));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &IDX(C, size, i, j));
			printcerr("", IDX(C, size, i, j), "\n");
#endif
		}
	}
}

void
mulCABtrans(size, C, A, B)
	size_t size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	size_t i, j, k;

	for (i = 0u; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &IDX(C, size, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(IDX(C, size, i, j)) = 0.0;
			Im(IDX(C, size, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(IDX(C, size, i, j), IDX(A, size, i, k), IDX(B, size, j, k));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &IDX(C, size, i, j));
			printcerr("", IDX(C, size, i, j), "\n");
#endif
		}
	}
}

void
transA(size, A)
	size_t size;
	complex_t *A;
{
	size_t i, j;
	complex_t t;
	
	for (i = 0u; i < size; i++)
		for (j = i; j < size; j++) {
			t = IDX(A, size, i, j);
			IDX(A, size, i, j) = IDX(A, size, j, i);
			IDX(A, size, j, i) = t;
		}
	
	return;
}

void
transAB(size, A, B)
	size_t size;
	complex_t *A;
	complex_t *B;
{
	size_t i, j;
	
	for (i = 0u; i < size; i++)
		for (j = 0u; j < size; j++)
			IDX(A, size, i, j) = IDX(B, size, j, i);
	
	return;
}

void
contrans(size, A)
	size_t size;
	complex_t *A;
{
	size_t i, j;
	complex_t t;
	
	for (i = 0u; i < size; i++)
		for (j = i; j < size; j++) {
			con2(t, IDX(A, size, i, j));
			con2(IDX(A, size, i, j), IDX(A, size, j, i));
			IDX(A, size, j, i) = t;
		}
	
	return;
}

void
contransAB(size, A, B)
	size_t size;
	complex_t *A;
	complex_t *B;
{
	size_t i, j;
	
	for (i = 0u; i < size; i++)
		for (j = 0u; j < size; j++)
			con2(IDX(A, size, i, j), IDX(B, size, j, i));
	
	return;
}
