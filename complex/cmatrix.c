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
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &IDX(C, i, 0u, size));
#endif
		for (j = 0u; j < size; j++) {
			Re(IDX(C, i, j, size)) = 0.0;
			Im(IDX(C, i, j, size)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(IDX(C, i, j, size), IDX(A, i, k, size), IDX(B, k, j, size));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &IDX(C, i, j, size));
			printcerr(IDX(C, i, j, size));
			fprintf(stderr, "\n");
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
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &IDX(C, i, 0u, size));
#endif
		for (j = 0u; j < size; j++) {
			Re(IDX(C, i, j, size)) = 0.0;
			Im(IDX(C, i, j, size)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(IDX(C, i, j, size), IDX(A, k, i, size), IDX(B, k, j, size));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &IDX(C, i, j, size));
			printcerr(IDX(C, i, j, size));
			fprintf(stderr, "\n");
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
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &IDX(C, i, 0u, size));
#endif
		for (j = 0u; j < size; j++) {
			Re(IDX(C, i, j, size)) = 0.0;
			Im(IDX(C, i, j, size)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(IDX(C, i, j, size), IDX(A, i, k, size), IDX(B, j, k, size));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &IDX(C, i, j, size));
			printcerr(IDX(C, i, j, size));
			fprintf(stderr, "\n");
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
			t = IDX(A, i, j, size);
			IDX(A, i, j, size) = IDX(A, j, i, size);
			IDX(A, j, i, size) = t;
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
			IDX(A, i, j, size) = IDX(B, j, i, size);
	
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
			con2(t, IDX(A, i, j, size));
			con2(IDX(A, i, j, size), IDX(A, j, i, size));
			IDX(A, j, i, size) = t;
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
			con2(IDX(A, i, j, size), IDX(B, j, i, size));
	
	return;
}
