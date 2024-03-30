#include "complex.h"
#include "cmatrix.h"
#include "clu.h"

/* INPUT: A - array of polongers to rows of a square matrix having dimension N
 *		Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *		The permutation matrix is not stored as a matrix, but in an longeger vector P of size N+1 
 *		containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
 *		where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S	
 */
int
LUPDecompose(N, A, P, Tol)
	long N;
	complex_t *A;
	long *P;
	double Tol;
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
			if ((absA = mag(MIJ(N, A, k, i))) > maxA) { 
				maxA = absA;
				imax = k;
			}
		
		if (maxA < Tol) return 0; /* failure, matrix is degenerate */
		
		if (imax != i) {
			/* pivoting P */
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			/* pivoting rows of A */
			for (j = 0l; j < N; j++) {
				t0 = MIJ(N, A, i, j);
				MIJ(N, A, i, j) = MIJ(N, A, imax, j);
				MIJ(N, A, imax, j) = t0;
			}
			
			/* counting pivots starting from N (for determinant) */
			P[N]++;
		}
		
		for (j = i + 1l; j < N; j++) {
			div2(MIJ(N, A, j, i), MIJ(N, A, i, i));
			
			for (k = i + 1l; k < N; k++) {
				mul3(t0, MIJ(N, A, j, i), MIJ(N, A, i, k));
				sub2(MIJ(N, A, j, k), t0);
			}
		}
	}
	
	return 1; /* decomposition done */
}

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
			if (P[i] == j) mk1(MIJ(N, Pm, i, j));
			else mk0(MIJ(N, Pm, i, j));
			
			if (i > j) {
				MIJ(N, L, i, j) = MIJ(N, A, i, j);
				mk0(MIJ(N, U, i, j));
			} else if (i == j) {
				mk1(MIJ(N, L, i, j));
				MIJ(N, U, i, j) = MIJ(N, A, i, j);
			} else {
				mk0(MIJ(N, L, i, j));
				MIJ(N, U, i, j) = MIJ(N, A, i, j);
			}
	}
	
	return;
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
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
			mul3(t0, MIJ(N, A, i, k), x[k]);
			sub2(x[i], t0);
		}
	}

	for (i = N - 1l; i >= 0l; i--) {
		for (k = i + 1l; k < N; k++) {
			mul3(t0, MIJ(N, A, i, k), x[k]);
			sub2(x[i], t0);
		}
		
		div2(x[i], MIJ(N, A, i, i));
	}
	
	return;
}

/* INPUT: A,P filled in LUPDecompose; N - dimension
 * OUTPUT: IA is the inverse of the initial matrix
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
			if (P[i] == j) mk1(MIJ(N, IA, i, j));
			else mk0(MIJ(N, IA, i, j));
			
			for (k = 0l; k < i; k++) {
				mul3(t0, MIJ(N, A, i, k), MIJ(N, IA, k, j));
				sub2(MIJ(N, IA, i, j), t0);
			}
		}

		for (i = N - 1l; i >= 0l; i--) {
			for (k = i + 1l; k < N; k++) {
				mul3(t0, MIJ(N, A, i, k), MIJ(N, IA, k, j));
				sub2(MIJ(N, IA, i, j), t0);
			}
			
			div2(MIJ(N, IA, i, j), MIJ(N, A, i, i));
		}
	}
	
	return;
}

/* INPUT: A,P filled in LUPDecompose; N - dimension. 
 * OUTPUT: Function returns the determinant of the initial matrix
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
		mul2(det, MIJ(N, A, i, i));
	
	if ((P[N] ^ N) & 1l) neg(det);
	return det;
}

