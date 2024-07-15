#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

/*
	INPUT:
		A - pointer to a square matrix, indexable by the MIJ macro
		N - dim(A)
		eps - lowest allowed degeneracy limit
	OUTPUT:
		A is changed in situ, after the procedure it contains a
		copy of both matrices L - E and U as A = (L - E) + U
		such that PA = LU. The vector P of size N + 1 contains a
		compressed form of the permutation matrix, with each element
		being the column in which the matrix has 1 on the
		corresponding row. The last element P[N] = S + N, where S is the
		number of row exchanges needed for the determinant computation,
		det(P) = (-1)^S. Also check LUPExtract().
	RETURNS:
		0 - failure (degeneracy limit reached)
		1 - A contains the LU decomposition in situ
*/
int
LUPDecompose(N, A, P, eps)
	long N;
	complex_t *A;
	long *P;
	double eps;
{
	long i, j, k, imax; 
	double maxA, absA;
	complex_t t0;
	
	for (i = 0l; i <= N; i++)
		P[i] = i; /* Unit permutation matrix, P[N] initialized with N */
	
	for (i = 0l; i < N; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++)
			if ((absA = magZ(MIJ(N, A, k, i))) > maxA) { 
				maxA = absA;
				imax = k;
			}
		
		if (maxA < eps) return 0; /* failure, matrix is degenerate */
		
		if (imax != i) {
			/* pivoting P */
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			/* pivoting rows of A */
			for (j = 0l; j < N; j++)
				xchgZ(MIJ(N, A, i, j), MIJ(N, A, imax, j));
			
			/* counting pivots starting from N (for determinant) */
			P[N]++;
		}
		
		for (j = i + 1l; j < N; j++) {
			divZ2(MIJ(N, A, j, i), MIJ(N, A, i, i));
			
			for (k = i + 1l; k < N; k++) {
				mulZ3(t0, MIJ(N, A, j, i), MIJ(N, A, i, k));
				subZ2(MIJ(N, A, j, k), t0);
			}
		}
	}
	
	return 1; /* decomposition done */
}

/*
	INPUT:
		A, P as filled by LUPDecompose
		N - dim(A)
	OUTPUT:
		Pm - the permutation matrix corresponding to P
		L - lower triangular
		U - upper triangular
*/
void
LUPExtract(N, Pm, L, U, A, P)
	long N;
	complex_t *Pm;
	complex_t *L;
	complex_t *U;
	complex_t *A;
	long *P;
{
	long i, j;
	
	for (i = 0l; i < N; i++)
		for (j = 0l; j < N; j++) {
			if (P[i] == j) mkZ1(MIJ(N, Pm, i, j));
			else mkZ0(MIJ(N, Pm, i, j));
			
			if (i > j) {
				MIJ(N, L, i, j) = MIJ(N, A, i, j);
				mkZ0(MIJ(N, U, i, j));
			} else if (i == j) {
				mkZ1(MIJ(N, L, i, j));
				MIJ(N, U, i, j) = MIJ(N, A, i, j);
			} else {
				mkZ0(MIJ(N, L, i, j));
				MIJ(N, U, i, j) = MIJ(N, A, i, j);
			}
	}
	
	return;
}

/*
	INPUT:
		A, P as filled by LUPDecompose
		N - dim(A)
		b - rhs vector
	OUTPUT:
		x - solution of of A*x=b
*/
void
LUPSolve(N, A, P, x, b)
	long N;
	complex_t *A;
	long *P;
	complex_t *x;
	complex_t *b;
{
	long i, k;
	complex_t t0;
	
	for (i = 0l; i < N; i++) {
		x[i] = b[P[i]];

		for (k = 0l; k < i; k++) {
			mulZ3(t0, MIJ(N, A, i, k), x[k]);
			subZ2(x[i], t0);
		}
	}

	for (i = N - 1l; i >= 0l; i--) {
		for (k = i + 1l; k < N; k++) {
			mulZ3(t0, MIJ(N, A, i, k), x[k]);
			subZ2(x[i], t0);
		}
		
		divZ2(x[i], MIJ(N, A, i, i));
	}
	
	return;
}

/*
	INPUT:
		A, P as filled by LUPDecompose
		N - dim(A)
		b - rhs vector
	OUTPUT:
		IA - A^-1
*/
void
LUPInvert(N, IA, A, P)
	long N;
	complex_t *IA;
	complex_t *A;
	long *P;
{
	long i, j, k;
	complex_t t0;
	
	for (j = 0l; j < N; j++) {
		for (i = 0l; i < N; i++) {
			if (P[i] == j) mkZ1(MIJ(N, IA, i, j));
			else mkZ0(MIJ(N, IA, i, j));
			
			for (k = 0l; k < i; k++) {
				mulZ3(t0, MIJ(N, A, i, k), MIJ(N, IA, k, j));
				subZ2(MIJ(N, IA, i, j), t0);
			}
		}

		for (i = N - 1l; i >= 0l; i--) {
			for (k = i + 1l; k < N; k++) {
				mulZ3(t0, MIJ(N, A, i, k), MIJ(N, IA, k, j));
				subZ2(MIJ(N, IA, i, j), t0);
			}
			
			divZ2(MIJ(N, IA, i, j), MIJ(N, A, i, i));
		}
	}
	
	return;
}

/*
	INPUT:
		A, P as filled by LUPDecompose
		N - dim(A)
	RETURNS:
		det(A)
*/
complex_t
LUPDeterminant(N, A, P)
	long N;
	complex_t *A;
	long *P;
{
	long i;
	complex_t det;
	
	det = MIJ(N, A, 0l, 0l);
	for (i = 1l; i < N; i++)
		mulZ2(det, MIJ(N, A, i, i));
	
	if ((P[N] ^ N) & 1l) negZ(det);
	
	return det;
}

