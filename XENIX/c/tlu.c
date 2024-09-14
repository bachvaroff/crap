#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

#define SIZE (4l)
#define SIZESQ (SIZE * SIZE)

int
main()
{
	complex_t A[SIZESQ], Adec[SIZESQ], Ainv[SIZESQ], AinvA[SIZESQ];
	complex_t Pm[SIZESQ], L[SIZESQ], U[SIZESQ], PmA[SIZESQ], LU[SIZESQ];
	complex_t Pminv[SIZESQ], PminvLU[SIZESQ];
	complex_t x[SIZE], b[SIZE], Ax[SIZE];
	complex_t Ainvb[SIZE];
	complex_t t0;
	long P[SIZE + 1l];
	int i, j;
	int r;
	
	mkZ(MIJ(SIZE, A, 0, 0), 4.0, 0.0);
	mkZ(MIJ(SIZE, A, 0, 1), 3.0, -1.0);
	mkZ(MIJ(SIZE, A, 0, 2), 3.0, 1.0);
	mkZ(MIJ(SIZE, A, 0, 3), -2.0, 3.5);
	mkZ(MIJ(SIZE, A, 1, 0), 6.0, 0.0);
	mkZ(MIJ(SIZE, A, 1, 1), 3.0, -2.0);
	mkZ(MIJ(SIZE, A, 1, 2), 3.0, 0.0);
	mkZ(MIJ(SIZE, A, 1, 3), 1.5, 2.5);
	mkZ(MIJ(SIZE, A, 2, 0), 3.0, 1.0);
	mkZ(MIJ(SIZE, A, 2, 1), 4.0, 0.0);
	mkZ(MIJ(SIZE, A, 2, 2), 3.0, -3.0);
	mkZ(MIJ(SIZE, A, 2, 3), 0.2, 0.3);
	mkZ(MIJ(SIZE, A, 3, 0), 4.1, 1.7);
	mkZ(MIJ(SIZE, A, 3, 1), 3.6, 1.9);
	mkZ(MIJ(SIZE, A, 3, 2), -2.3, 0.25);
	mkZ(MIJ(SIZE, A, 3, 3), -5.0, -1.15);
	
	mkZ(b[0], 1.0, 0.0);
	mkZ(b[1], 0.0, 1.0);
	mkZ(b[2], 0.5, 0.5);
	mkZ(b[3], 2.2, 3.4);
	
	cpAB(SIZE, Adec, A);
	
	printm(SIZE, A, "A =\n", 0, 0);
	r = LUPDecompose(SIZE, Adec, P, 0.000000000001);
	assert(r);
	printm(SIZE, Adec, "LU(A) =\n", 0, 0);
	for (i = 0; i <= SIZE; i++)
		printf("P[%d] = %ld\n", i, P[i]);
	
	LUPExtract(SIZE, Pm, L, U, Adec, P);
	printm(SIZE, Pm, "Pm =\n", 0, 0);
	printm(SIZE, L, "L =\n", 0, 0);
	printm(SIZE, U, "U =\n", 0, 0);
	
	mulCAB(SIZE, PmA, Pm, A);
	mulCAB(SIZE, LU, L, U);
	printm(SIZE, PmA, "Pm*A =\n", 0, 0);
	printm(SIZE, LU, "L*U =\n", 0, 0);
	
	cpActB(SIZE, Pminv, Pm);
	mulCAB(SIZE, PminvLU, Pminv, LU);
	printm(SIZE, PminvLU, "inv(Pm)*L*U =\n", 0, 0);
	
	t0 = LUPDeterminant(SIZE, Adec, P);
	printf("det(A) = ");
	printZ("", t0, "\n");
	
	LUPInvert(SIZE, Ainv, Adec, P);
	printm(SIZE, Ainv, "inv(A) =\n", 0, 0);
	
	mulCAB(SIZE, AinvA, Ainv, A);
	printm(SIZE, AinvA, "inv(A)*A =\n", 0, 0);
	
	LUPSolve(SIZE, Adec, P, x, b);
	printv(SIZE, x, "x =\n", 0);
	
	mulYvAXv(SIZE, Ax, A, x);
	printv(SIZE, Ax, "A*x =\n", 0);
	printv(SIZE, b, "b =\n", 0);
	
	mulYvAXv(SIZE, Ainvb, Ainv, b);
	printv(SIZE, Ainvb, "inv(A)*b =\n", 0);
	
	return 0;
}
