#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"
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
		fflush(stderr);
#endif
		for (j = 0l; j < size; j++) {
			mkZ0(MIJ(size, C, i, j));
			for (k = 0l; k < size; k++)
				maddZ(MIJ(size, C, i, j),
					MIJ(size, A, i, k),
					MIJ(size, B, k, j));
#ifdef _DEBUG_CI_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_
		if (j & 0x7l) fputc((int)'_', stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
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
			mkZ0(MIJ(size, C, i, j));
			for (k = 0l; k < size; k++)
				maddZ(MIJ(size, C, i, j),
					MIJ(size, A, k, i),
					MIJ(size, B, k, j));
#ifdef _DEBUG_CI_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_
		if (j & 0x7l) fputc((int)'_', stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
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
			mkZ0(MIJ(size, C, i, j));
			for (k = 0l; k < size; k++)
				maddZ(MIJ(size, C, i, j),
					MIJ(size, A, i, k),
					MIJ(size, B, j, k));
#ifdef _DEBUG_CI_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_
		if (j & 0x7l) fputc((int)'_', stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
	}
}

void
mulCAcontransB(size, C, A, B)
	long size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	long i, j, k;
	complex_t t0;
	
	for (i = 0l; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%16lx % 11ld %16lx\n",
			(unsigned long)C, i,
			(unsigned long)&MIJ(size, C, i, 0l));
#endif
		for (j = 0l; j < size; j++) {
			mkZ0(MIJ(size, C, i, j));
			for (k = 0l; k < size; k++) {
				conZ2(t0, MIJ(size, A, k, i));
				maddZ(MIJ(size, C, i, j),
					t0, MIJ(size, B, k, j));
			}
#ifdef _DEBUG_CI_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_
		if (j & 0x7l) fputc((int)'_', stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
	}
}

void
mulCABcontrans(size, C, A, B)
	long size;
	complex_t *C;
	complex_t *A;
	complex_t *B;
{
	long i, j, k;
	complex_t t0;
	
	for (i = 0l; i < size; i++) {
#ifdef _DEBUG_CI_
		fprintf(stderr, "%16lx % 11ld %16lx\n",
			(unsigned long)C, i,
			(unsigned long)&MIJ(size, C, i, 0l));
#endif
		for (j = 0l; j < size; j++) {
			mkZ0(MIJ(size, C, i, j));
			for (k = 0l; k < size; k++) {
				conZ2(t0, MIJ(size, B, j, k));
				maddZ(MIJ(size, C, i, j),
					MIJ(size, A, i, k), t0);
			}
#ifdef _DEBUG_CI_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_
		if (j & 0x7l) fputc((int)'_', stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
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
			conZ2(t, MIJ(size, A, i, j));
			conZ2(MIJ(size, A, i, j), MIJ(size, A, j, i));
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
			conZ2(MIJ(size, A, i, j), MIJ(size, B, j, i));
	
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
		mkZ0(y[i]);
		for (j = 0l; j < size; j++)
			maddZ(y[i], MIJ(size, A, i, j), x[j]);
	}
	
	return;
}

