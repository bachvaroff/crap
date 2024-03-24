#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

#define SIZE 20l
#define SIZESQ (SIZE * SIZE)

complex_t *a, *b, *dec, *abi, *biai;
long P[SIZE + 1l];

complex_t *marr;

int
main()
{
	long i, j;
	complex_t t0;
	double t1;
	double maxij;
	long maxi, maxj;
	
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
			if (i == j) expiphi(MIJ(SIZE, b, i, j),
				M_PI / 3.0 + 2.0 * M_PI * (double)i / (double)SIZE);
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
			sub3(t0, MIJ(SIZE, abi, i, j), MIJ(SIZE, biai, i, j));
			t1 = log(mag(t0)) / M_LN10;
			if (i || j) {
				if (maxij < t1) {
					maxij = t1;
					maxi = i;
					maxj = j;
				}
			} else {
				maxij = t1;
				maxi = i;
				maxj = j;
			}
		}
	printf("%ld | %ld | %lf\n", maxi, maxj, maxij);
	printf("abi[%ld][%ld] = ", maxi, maxj);
	printc("", MIJ(SIZE, abi, maxi, maxj), " | ");
	printf("biai[%ld][%ld] = ", maxi, maxj);
	printc("", MIJ(SIZE, biai, maxi, maxj), "\n");
	
	printf("end\n");

	free((void *)marr);
	
bad0:
	return 0;
}
