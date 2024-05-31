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
#ifdef _DEBUG_CI_D_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_D_
		j &= 0x7l;
		if (j) fputc((int)'0' + (int)j, stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
	}
}

void
mulCtransAB(size, C, A, B)
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
#ifdef _DEBUG_CI_D_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_D_
		j &= 0x7l;
		if (j) fputc((int)'0' + (int)j, stderr);
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
					MIJ(size, A, i, k),
					MIJ(size, B, j, k));
#ifdef _DEBUG_CI_D_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_D_
		j &= 0x7l;
		if (j) fputc((int)'0' + (int)j, stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
	}
}

void
mulCconjtransAB(size, C, A, B)
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
				conjZ2(t0, MIJ(size, A, k, i));
				maddZ(MIJ(size, C, i, j),
					t0, MIJ(size, B, k, j));
			}
#ifdef _DEBUG_CI_D_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_D_
		j &= 0x7l;
		if (j) fputc((int)'0' + (int)j, stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
	}
}

void
mulCAconjtransB(size, C, A, B)
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
				conjZ2(t0, MIJ(size, B, j, k));
				maddZ(MIJ(size, C, i, j),
					MIJ(size, A, i, k), t0);
			}
#ifdef _DEBUG_CI_D_
			if ((j & 0x7l) == 0x7l) {
				fputc((int)'.', stderr);
				fflush(stderr);
			}
#endif
		}
#ifdef _DEBUG_CI_D_
		j &= 0x7l;
		if (j) fputc((int)'0' + (int)j, stderr);
		fputc((int)'\n', stderr);
		fflush(stderr);
#endif
	}
}

void
AconjA(size, A)
	long size; 
	complex_t *A;
{
	long i, j;
	
	for (i = 0l; i < size; i++)
		for (j = i; j < size; j++)
                        conjZ(MIJ(size, A, i, j));
	
	return;
}

void  
AconjB(size, A, B)
	long size;
	complex_t *A;
	complex_t *B;
{
	long i, j;
	
	for (i = 0l; i < size; i++)
		for (j = i; j < size; j++)
			conjZ2(MIJ(size, A, i, j), MIJ(size, B, i, j));
        
        return;
}

void
AtransA(size, A)
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
AtransB(size, A, B)
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
AconjtransA(size, A)
	long size;
	complex_t *A;
{
	long i, j;
	complex_t t;
	
	for (i = 0l; i < size; i++)
		for (j = i; j < size; j++) {
			conjZ2(t, MIJ(size, A, i, j));
			conjZ2(MIJ(size, A, i, j), MIJ(size, A, j, i));
			MIJ(size, A, j, i) = t;
		}
	
	return;
}

void
AconjtransB(size, A, B)
	long size;
	complex_t *A;
	complex_t *B;
{
	long i, j;
	
	for (i = 0l; i < size; i++)
		for (j = 0l; j < size; j++)
			conjZ2(MIJ(size, A, i, j), MIJ(size, B, j, i));
	
	return;
}

void
mulyvAxv(size, yv, A, xv)
	long size;
	complex_t *yv;
	complex_t *A;
	complex_t *xv;
{
	long i, j;
	
	for (i = 0l; i < size; i++) {
		mkZ0(yv[i]);
		for (j = 0l; j < size; j++)
			maddZ(yv[i], MIJ(size, A, i, j), xv[j]);
	}
	
	return;
}

void
mulycovxcovA(size, ycv, xcv, A)
	long size;
	complex_t *ycv;
	complex_t *xcv;
	complex_t *A;
{
	long i, j;
	
	for (j = 0l; j < size; j++) {
		mkZ0(ycv[j]);
		for (i = 0l; i < size; i++)
			maddZ(ycv[j], xcv[i], MIJ(size, A, i, j));
	}
	
	return;
}
