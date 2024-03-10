#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

int
main()
{
	complex_t z0, z1, u;

	mkC(z0, 1.0, 0.5);
	z1 = z0;
	
	CosZ(z0);
	SinZ(z1);
	printc("", z0, "\n");
	printc("", z1, "\n");
	mul2(z0, z0);
	mul2(z1, z1);
	printc("", z0, "\n");
	printc("", z1, "\n");
	add3(u, z0, z1);
	printc("", u, "\n");
	
	LogZ(z0);
	printc("", z0, "\n");
	
	expZ(z0);
	printc("", z0, "\n");
	
	LogZ(z1);
	printc("", z1, "\n");
	
	expZ(z1);
	printc("", z1, "\n");

	return 0;
}
