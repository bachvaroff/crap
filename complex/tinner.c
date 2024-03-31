#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

int
main()
{
	complex_t z1, z2, z3, z4, z5, t0, t1, t2;
	
	Re(z1) = 1.0;
	Im(z1) = 0.5;
	
	Re(z2) = 2.5;
	Im(z2) = -1.5;
	
	innerZ3(t0, z1, z2);
	printZ("", t0, "\n");
	
	innerZ3(t0, z2, z1);
	printZ("", t0, "\n");
	
	conZ(t0);
	printZ("", t0, "\n");
	
	Re(z3) = 2.0;
	Im(z3) = -sqrt(2.0) / 2.0;

	Re(z4) = sqrt(3.0);
	Im(z4) = -sqrt(2.0) / 2.0;
	
	Re(z5) = Im(z5) = 1.0;
	
	mulZ3(t0, z1, z3);
	maddZ(t0, z2, z4);
	innerZ3(t1, t0, z5);
	printZ("", t1, "\n");
	
	innerZ3(t0, z3, z5);
	mulZ2(t0, z1);
	innerZ3(t1, z4, z5);
	maddZ(t0, t1, z2);
	printZ("", t0, "\n");
	
	mulZ3(t0, z1, z3);
	maddZ(t0, z2, z4);
	innerZ3(t1, z5, t0);
	printZ("", t1, "\n");

	innerZ3(t0, z5, z3);
	conZ2(t1, z1);
	mulZ2(t0, t1);
	innerZ3(t1, z5, z4);
	conZ2(t2, z2);
	maddZ(t0, t1, t2);
	printZ("", t0, "\n");
	
	return 0;
}
