#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

int
main()
{
	complex_t A[9], decA[9], IA[9], C[9];
	complex_t x[3], b[3], y[3];
	complex_t t0;
	long P[4];
	int i, j;
	int r;
	
	mkC(MIJ(3, A, 0, 0), 4.0, 0.0);
	mkC(MIJ(3, A, 0, 1), 3.0, -1.0);
	mkC(MIJ(3, A, 0, 2), 3.0, 1.0);
	mkC(MIJ(3, A, 1, 0), 6.0, 0.0);
	mkC(MIJ(3, A, 1, 1), 3.0, -2.0);
	mkC(MIJ(3, A, 1, 2), 3.0, 0.0);
	mkC(MIJ(3, A, 2, 0), 3.0, 1.0);
	mkC(MIJ(3, A, 2, 1), 4.0, 0.0);
	mkC(MIJ(3, A, 2, 2), 3.0, -3.0);
	
	mkC(b[0], 1.0, 0.0);
	mkC(b[1], 0.0, 1.0);
	mkC(b[2], 0.5, 0.5);
	
	copyAB(3, decA, A);
	
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("A[%d][%d] = ", i, j);
			printc("", MIJ(3, A, i, j), "\n");
		}
	r = LUPDecompose(3, decA, P, 0.000000001);
	assert(r);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("\tdecA[%d][%d] = ", i, j);
			printc("", MIJ(3, decA, i, j), "\n");
		}
	for (i = 0; i <= 3; i++)
		printf("\tP[%d] = %ld\n", i, P[i]);
	
	t0 = LUPDeterminant(3, decA, P);
	printf("det(A) = ");
	printc("", t0, "\n");
	
	LUPInvert(3, IA, decA, P);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("IA[%d][%d] = ", i, j);
			printc("", MIJ(3, IA, i, j), "\n");
		}
	
	mulCAB(3, C, A, IA);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("C[%d][%d] = ", i, j);
			printc("", MIJ(3, C, i, j), "\n");
		}
	
	LUPSolve(3, decA, P, x, b);
	for (i = 0; i < 3; i++) {
		printf("x[%d] =", i);
		printc("", x[i], "\n");
	}
	
	mulyAx(3, y, A, x);
	for (i = 0; i < 3; i++) {
		printf("y[%d] =", i);
		printc("", y[i], " | ");
		printf("b[%d] =", i);
		printc("", b[i], "\n");
	}
	
	return 0;
}
