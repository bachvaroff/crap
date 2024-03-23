#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

#define SIZE 400l
#define SIZESQ (SIZE * SIZE)

complex_t *a, *b, *dec, *abi, *biai;
long P[SIZE + 1l];

complex_t *marr;

int
main()
{
	long i, j;
	complex_t t0;
	
	printf("init\n");
	
	marr = (complex_t *)calloc(5l * SIZESQ, sizeof (complex_t));
	if (!marr) {
		fprintf(stderr, "cannot alloc marr\n");
		goto bad0;
	} else printf("marr %16lx % 10ld\n", (unsigned long)marr, 5l * SIZESQ * sizeof (complex_t));

	a = marr + 0l * SIZESQ;
	b = marr + 1l * SIZESQ;
	dec = marr + 2l * SIZESQ;
	abi = marr + 3l * SIZESQ;
	biai = marr + 4l * SIZESQ;
	printf("a %16lx | b %16lx | dec %16lx | abi %16lx | biai %16lx\n",
		(unsigned long)a, (unsigned long)b, (unsigned long)dec,
		(unsigned long)abi, (unsigned long)biai);

	for (i = 0l; i < SIZE; i++)
		for (j = 0l; j < SIZE; j++) {
			if (i == j) expiphi(MIJ(SIZE, a, i, j),
				2.0 * M_PI * (double)i / (double)SIZE);
			else mkC(MIJ(SIZE, a, i, j),
				cos(M_PI * (double)j / (double)SIZE),
				sin(M_PI * (double)i / (double)SIZE));
		}
	
	for (i = 0l; i < SIZE; i++)
		for (j = 0l; j < SIZE; j++) {
			if (i == j) mkC(MIJ(SIZE, b, i, j), 1.0, 0.0);
			else mkC(MIJ(SIZE, b, i, j),
				cos(M_PI * (double)i / (double)SIZE),
				sin(M_PI * (double)j / (double)SIZE));
		}
	
	printf("begin\n");
	
	mulCAB(SIZE, abi, a, b);
	printf("ab\n");
	
	copyAB(SIZE, dec, abi);
	assert(LUPDecompose(SIZE, dec, P, 0.000000001));
	LUPInvert(SIZE, abi, dec, P);
	printf("abi\n");
	
	copyAB(SIZE, dec, a);
	assert(LUPDecompose(SIZE, dec, P, 0.000000001));
	LUPInvert(SIZE, a, dec, P);
	printf("ai\n");
	
	copyAB(SIZE, dec, b);
	assert(LUPDecompose(SIZE, dec, P, 0.000000001));
	LUPInvert(SIZE, b, dec, P);
	printf("bi\n");
	
	mulCAB(SIZE, biai, b, a);
	printf("biai\n");
	
	for (i = 0l; i < SIZE; i++)
		for (j = 0l; j < SIZE; j++) {
/*
			printf("abi[%ld][%ld] = ", i, j);
			printc("", MIJ(SIZE, abi, i, j), " | ");
			printf("biai[%ld][%ld] = ", i, j);
			printc("", MIJ(SIZE, biai, i, j), " | ");
*/
			sub3(t0, MIJ(SIZE, abi, i, j), MIJ(SIZE, biai, i, j));
			printc("", t0, " | ");
			printf("%lf\n", log(mag(t0)) / M_LN10);
		}
	
	printf("end\n");

	free((void *)marr);
	
bad0:
	return 0;
}
