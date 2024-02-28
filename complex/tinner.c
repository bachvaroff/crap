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
	
	inner3(t0, z1, z2);
	printc("", t0, "\n");
	
	inner3(t0, z2, z1);
	printc("", t0, "\n");
	
	con(t0);
	printc("", t0, "\n");
	
	Re(z3) = 2.0;
	Im(z3) = -sqrt(2.0) / 2.0;

	Re(z4) = sqrt(3.0);
	Im(z4) = -sqrt(2.0) / 2.0;
	
	Re(z5) = Im(z5) = 1.0;
	
	mul3(t0, z1, z3);
	madd(t0, z2, z4);
	inner3(t1, t0, z5);
	printc("", t1, "\n");
	
	inner3(t0, z3, z5);
	mul2(t0, z1);
	inner3(t1, z4, z5);
	madd(t0, t1, z2);
	printc("", t0, "\n");
	
	mul3(t0, z1, z3);
	madd(t0, z2, z4);
	inner3(t1, z5, t0);
	printc("", t1, "\n");

	inner3(t0, z5, z3);
	con2(t1, z1);
	mul2(t0, t1);
	inner3(t1, z5, z4);
	con2(t2, z2);
	madd(t0, t1, t2);
	printc("", t0, "\n");
	
	return 0;
}
