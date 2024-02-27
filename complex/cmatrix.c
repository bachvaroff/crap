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
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &MIJ(C, size, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(MIJ(C, size, i, j)) = 0.0;
			Im(MIJ(C, size, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(MIJ(C, size, i, j), MIJ(A, size, i, k), MIJ(B, size, k, j));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &MIJ(C, size, i, j));
			printcerr("", MIJ(C, size, i, j), "\n");
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
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &MIJ(C, size, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(MIJ(C, size, i, j)) = 0.0;
			Im(MIJ(C, size, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(MIJ(C, size, i, j), MIJ(A, size, k, i), MIJ(B, size, k, j));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &MIJ(C, size, i, j));
			printcerr("", MIJ(C, size, i, j), "\n");
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
		fprintf(stderr, "%08x % 11u %08x\n", C, i, &MIJ(C, size, i, 0u));
#endif
		for (j = 0u; j < size; j++) {
			Re(MIJ(C, size, i, j)) = 0.0;
			Im(MIJ(C, size, i, j)) = 0.0;
			for (k = 0u; k < size; k++)
				madd(MIJ(C, size, i, j), MIJ(A, size, i, k), MIJ(B, size, j, k));
#ifdef _DEBUG_CIJ_
			fprintf(stderr, "\t%08x % 11u% 11u %08x ", C, i, j, &MIJ(C, size, i, j));
			printcerr("", MIJ(C, size, i, j), "\n");
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
			t = MIJ(A, size, i, j);
			MIJ(A, size, i, j) = MIJ(A, size, j, i);
			MIJ(A, size, j, i) = t;
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
			MIJ(A, size, i, j) = MIJ(B, size, j, i);
	
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
			con2(t, MIJ(A, size, i, j));
			con2(MIJ(A, size, i, j), MIJ(A, size, j, i));
			MIJ(A, size, j, i) = t;
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
			con2(MIJ(A, size, i, j), MIJ(B, size, j, i));
	
	return;
}
