#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

int
main()
{
	complex_t z0, z1, u;

	mkZ(z0, 1.0, 0.5);
	z1 = z0;
	
	CosZ(z0);
	SinZ(z1);
	printZ("", z0, "\n");
	printZ("", z1, "\n");
	mulZ2(z0, z0);
	mulZ2(z1, z1);
	printZ("", z0, "\n");
	printZ("", z1, "\n");
	addZ3(u, z0, z1);
	printZ("", u, "\n");
	
	LogZ(z0);
	printZ("", z0, "\n");
	
	expZ(z0);
	printZ("", z0, "\n");
	
	LogZ(z1);
	printZ("", z1, "\n");
	
	expZ(z1);
	printZ("", z1, "\n");
	
	TanZ(z0);
	printZ("", z0, "\n");
	
	TanZ(z1);
	printZ("", z1, "\n");
	
	return 0;
}
