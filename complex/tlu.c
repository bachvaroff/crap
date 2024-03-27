#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

int
main()
{
	complex_t A[9], Adec[9], Ainv[9], AAinv[9];
	complex_t Pm[9], L[9], U[9], PA[9], LU[9];
	complex_t x[3], b[3], Ax[3];
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
	
	copyAB(3, Adec, A);
	
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("A[%d][%d] = ", i, j);
			printc("", MIJ(3, A, i, j), "\n");
		}
	r = LUPDecompose(3, Adec, P, 0.000000001);
	assert(r);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("\tAdec[%d][%d] = ", i, j);
			printc("", MIJ(3, Adec, i, j), "\n");
		}
	for (i = 0; i <= 3; i++)
		printf("\tP[%d] = %ld\n", i, P[i]);
	
	LUPExtract(3, Pm, L, U, Adec, P);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("Pm[%d][%d] = ", i, j);
			printc("", MIJ(3, Pm, i, j), " | ");
			printf("L[%d][%d] = ", i, j);
			printc("", MIJ(3, L, i, j), " | ");
			printf("U[%d][%d] = ", i, j);
			printc("", MIJ(3, U, i, j), "\n");
		}
	mulCAB(3, PA, Pm, A);
	mulCAB(3, LU, L, U);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("P*A[%d][%d] = ", i, j);
			printc("", MIJ(3, PA, i, j), " | ");
			printf("L*U[%d][%d] = ", i, j);
			printc("", MIJ(3, LU, i, j), "\n");
		}
	
	t0 = LUPDeterminant(3, Adec, P);
	printf("detA = ");
	printc("", t0, "\n");
	
	LUPInvert(3, Ainv, Adec, P);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("Ainv[%d][%d] = ", i, j);
			printc("", MIJ(3, Ainv, i, j), "\n");
		}
	
	mulCAB(3, AAinv, A, Ainv);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			printf("A*Ainv[%d][%d] = ", i, j);
			printc("", MIJ(3, AAinv, i, j), "\n");
		}
	
	LUPSolve(3, Adec, P, x, b);
	for (i = 0; i < 3; i++) {
		printf("x[%d] =", i);
		printc("", x[i], "\n");
	}
	
	mulyAx(3, Ax, A, x);
	for (i = 0; i < 3; i++) {
		printf("A*x[%d] =", i);
		printc("", Ax[i], " | ");
		printf("b[%d] =", i);
		printc("", b[i], "\n");
	}
	
	return 0;
}
