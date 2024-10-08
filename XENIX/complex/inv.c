#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

#define SIZE (400l)
#define SIZESQ (SIZE * SIZE)

complex_t *a, *b, *dec, *abi, *biai;
long P[SIZE + 1l];

complex_t *marr;

int
main()
{
	long i, j;
	complex_t t0;
	double mi, stddev;
	
	printf("init(N = %ld)\n", SIZE);
	
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
			if (i == j) expiphiZ(MIJ(SIZE, a, i, j),
				2.0 * M_PI * (double)i / (double)SIZE);
			else mkZ(MIJ(SIZE, a, i, j),
				0.0,
				sin(M_PI * (double)i / (double)SIZE));
		}
	
	for (i = 0l; i < SIZE; i++)
		for (j = 0l; j < SIZE; j++) {
			if (i == j) expiphiZ(MIJ(SIZE, b, i, j),
				2.0 * M_PI * (double)i / (double)SIZE);
			else mkZ(MIJ(SIZE, b, i, j),
				cos(M_PI * (double)j / (double)SIZE),
				0.0);
		}
	
	printf("begin\n");
	
	mulCAB(SIZE, abi, a, b);
	cpAB(SIZE, dec, abi);
	assert(LUPDecompose(SIZE, dec, P, 0.000000000001));
	t0 = LUPDeterminant(SIZE, dec, P);
	printZ("det(ab) = ", t0, "\n");
	invZ(t0);
	printZ("1/det(ab) = ", t0, "\n**** ab ****\n");
	
	LUPInvert(SIZE, abi, dec, P);
	cpAB(SIZE, dec, abi);
	assert(LUPDecompose(SIZE, dec, P, 0.000000000001));
	t0 = LUPDeterminant(SIZE, dec, P);
	printZ("det(inv(ab)) = ", t0, "\n**** inv(ab) ****\n");
	
	cpAB(SIZE, dec, a);
	assert(LUPDecompose(SIZE, dec, P, 0.000000000001));
	t0 = LUPDeterminant(SIZE, dec, P);
	printZ("det(a) = ", t0, "\n");
	invZ(t0);
	printZ("1/det(a) = ", t0, "\n");
	LUPInvert(SIZE, a, dec, P);
	cpAB(SIZE, dec, a);
	assert(LUPDecompose(SIZE, dec, P, 0.000000000001));
	t0 = LUPDeterminant(SIZE, dec, P);
	printZ("det(inv(a)) = ", t0, "\n**** inv(a) ****\n");
	
	cpAB(SIZE, dec, b);
	assert(LUPDecompose(SIZE, dec, P, 0.000000000001));
	t0 = LUPDeterminant(SIZE, dec, P);
	printZ("det(b) = ", t0, "\n");
	invZ(t0);
	printZ("1/det(b) = ", t0, "\n");
	LUPInvert(SIZE, b, dec, P);
	cpAB(SIZE, dec, b);
	assert(LUPDecompose(SIZE, dec, P, 0.000000000001));
	t0 = LUPDeterminant(SIZE, dec, P);
	printZ("det(inv(b)) = ", t0, "\n**** inv(b) ****\n");
	
	mulCAB(SIZE, biai, b, a);
	cpAB(SIZE, dec, biai);
	assert(LUPDecompose(SIZE, dec, P, 0.000000000001));
	t0 = LUPDeterminant(SIZE, dec, P);
	printZ("det(inv(b)inv(a)) = ", t0, "\n**** inv(b)inv(a) ****\n");
	
	mi = 0.0;
	for (i = 0l; i < SIZE; i++)
		for (j = 0l; j < SIZE; j++) {
			subZ3(t0, MIJ(SIZE, abi, i, j), MIJ(SIZE, biai, i, j));
			mi += magZ(t0);
		}
	mi /= (double)SIZESQ;
	
	stddev = 0.0;
	for (i = 0l; i < SIZE; i++)
		for (j = 0l; j < SIZE; j++) {
			subZ3(t0, MIJ(SIZE, abi, i, j), MIJ(SIZE, biai, i, j));
			stddev += (magZ(t0) - mi) * (magZ(t0) - mi);
		}
	stddev = sqrt(stddev / (double)SIZESQ);
	
	printf("mi = %.16lf | stddev = %.16lf\n", mi, stddev);
		
	printf("end\n");

	free((void *)marr);
	
bad0:
	return 0;
}
