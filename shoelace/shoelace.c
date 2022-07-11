#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#define PERM4 24

typedef struct point {
	float X, Y;
} point_t;

static const int perm[PERM4][4] = {
	{ 0, 1, 2, 3 },
	{ 0, 1, 3, 2 },
	{ 0, 2, 1, 3 },
	{ 0, 2, 3, 1 },
	{ 0, 3, 1, 2 },
	{ 0, 3, 2, 1 },
	{ 1, 0, 2, 3 },
	{ 1, 0, 3, 2 },
	{ 1, 2, 0, 3 },
	{ 1, 2, 3, 0 },
	{ 1, 3, 0, 2 },
	{ 1, 3, 2, 0 },
	{ 2, 0, 1, 3 },
	{ 2, 0, 3, 1 },
	{ 2, 1, 0, 3 },
	{ 2, 1, 3, 0 },
	{ 2, 3, 0, 1 },
	{ 2, 3, 1, 0 },
	{ 3, 0, 1, 2 },
	{ 3, 0, 2, 1 },
	{ 3, 1, 0, 2 },
	{ 3, 1, 2, 0 },
	{ 3, 2, 0, 1 },
	{ 3, 2, 1, 0 }
};

float shoelace(point_t A, point_t B, point_t C, point_t D) {
	float x1 = A.X, y1 = A.Y;
	float x2 = B.X, y2 = B.Y;
	float x3 = C.X, y3 = C.Y;
	float x4 = D.X, y4 = D.Y;
	
	return 0.5f * (
			(x1 * y2 + x2 * y3 + x3 * y4 + x4 * y1) - 
			(x2 * y1 + x3 * y2 + x4 * y3 + x1 * y4)
	);
}

int cross(point_t A, point_t B, point_t C) {
	float x1 = A.X, y1 = A.Y;
	float x2 = B.X, y2 = B.Y;
	float x3 = C.X, y3 = C.Y;
	
	return (y3 - y1) * (x2 - x1) > (y2 - y1) * (x3 - x1);
}

int intersect(point_t A, point_t B, point_t C, point_t D) {
	    return (cross(A, C, D) != cross(B, C, D)) && (cross(A, B, C) != cross(A, B, D));
}

int main(void) {
	int j;
	point_t temp[4];
	
	scanf("%f%f%f%f%f%f%f%f",
			&temp[0].X, &temp[0].Y,
			&temp[1].X, &temp[1].Y,
			&temp[2].X, &temp[2].Y,
			&temp[3].X, &temp[3].Y
	);
	
	for (j = 0; j < PERM4; j++)
		printf("%d %f\n",
				intersect(temp[perm[j][0]], temp[perm[j][1]], temp[perm[j][2]], temp[perm[j][3]]),
				shoelace(temp[perm[j][0]], temp[perm[j][1]], temp[perm[j][2]], temp[perm[j][3]])
		);
	
	return 0;
}

