#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"
#include "cmatrix.h"

void
printm(size, A, label, pr, pc)
	long size;
	complex_t *A;
	char *label;
	int pr;
	int pc;
{
	long i, j;
	
	if (label) printf("%s", label);
	for (i = 0l; i < size; i++) {
		fputc('\t', stdout);
		if (pr) printf("[% 11ld] | ", i);
		for (j = 0l; j < size; j++) {
			if (pc) printf("[% 11ld]", j);
			printZ("", MIJ(size, A, i, j), "");
		}
		fputc('\n', stdout);
	}
	
	return;
}

void
printv(size, xv, label, pr)
	long size;
	complex_t *xv;
	char *label;
	int pr;
{
	long i;
	
	if (label) printf("%s", label);
	for (i = 0l; i < size; i++) {
		fputc('\t', stdout);
		if (pr) printf("[% 11ld] | ", i);
		printZ("", xv[i], "\n");
	}
	
	return;
}

void
printcv(size, xcv, label, pc)
	long size;
	complex_t *xcv;
	char *label;
	int pc;
{
	long i;
	
	if (label) printf("%s", label);
	fputc('\t', stdout);
	for (i = 0l; i < size; i++) {
		if (pc) printf("[% 11ld]", i);
		printZ("", xcv[i], "");
	}
	fputc('\n', stdout);
	
	return;
}

void
isAcA(size, A)
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
cpAcB(size, A, B)
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
isAtA(size, A)
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
cpAtB(size, A, B)
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
isActA(size, A)
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
cpActB(size, A, B)
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
mulCtAB(size, C, A, B)
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
mulCctAB(size, C, A, B)
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
mulCAtB(size, C, A, B)
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
mulCActB(size, C, A, B)
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
mulYvAXv(size, yv, A, xv)
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
mulYcvXcvA(size, ycv, xcv, A)
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

