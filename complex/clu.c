#include "cmatrix.h"
#include "clu.h"

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *		Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *		The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
 *		containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
 *		where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S	
 */
int
LUPDecompose(N, A, P, Tol)
	int N;
	complex_t *A;
	int *P;
	double Tol;
{
	int i, j, k, imax; 
	double maxA, absA;
	complex_t t0;
	
	for (i = 0; i <= N; i++)
		P[i] = i; /* Unit permutation matrix, P[N] initialized with N */
	
	for (i = 0; i < N; i++) {
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
			for (j = 0; j < N; j++) {
				t0 = MIJ(N, A, i, j);
				MIJ(N, A, i, j) = MIJ(N, A, imax, j);
				MIJ(N, A, imax, j) = t0;
			}
			
			/* counting pivots starting from N (for determinant) */
			P[N]++;
		}
		
		for (j = i + 1; j < N; j++) {
			div2(MIJ(N, A, j, i), MIJ(N, A, i, i));
			
			for (k = i + 1; k < N; k++) {
				mul3(t0, MIJ(N, A, j, i), MIJ(N, A, i, k));
				sub2(MIJ(N, A, j, k), t0);
			}
		}
	}
	
	return 1; /* decomposition done */
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void
LUPSolve(N, A, P, x, b)
	int N;
	complex_t *A;
	int *P;
	complex_t *x;
	complex_t *b;
{
	int i, k;
	complex_t t0;
	
	for (i = 0; i < N; i++) {
		x[i] = b[P[i]];

		for (k = 0; k < i; k++) {
			mul3(t0, MIJ(N, A, i, k), x[k]);
			sub2(x[i], t0);
		}
	}

	for (i = N - 1; i >= 0; i--) {
		for (k = i + 1; k < N; k++) {
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
	int N;
	complex_t *IA;
	complex_t *A;
	int *P;
{
	int i, j, k;
	complex_t t0;
	
	for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {
			if (P[i] == j) mk1(MIJ(N, IA, i, j));
			else mk0(MIJ(N, IA, i, j));
			
			for (k = 0; k < i; k++) {
				mul3(t0, MIJ(N, A, i, k), MIJ(N, IA, k, j));
				sub2(MIJ(N, IA, i, j), t0);
			}
		}

		for (i = N - 1; i >= 0; i--) {
			for (k = i + 1; k < N; k++) {
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
	int N;
	complex_t *A;
	int *P;
{
	int i;
	complex_t det;
	
	det = MIJ(N, A, 0, 0);

	for (i = 1; i < N; i++)
		mul2(det, MIJ(N, A, i, i));
	
	if ((P[N] - N) % 2) scale2(det, -1.0);
	return det;
}

