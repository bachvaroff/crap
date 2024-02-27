#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

int
main()
{
	complex_t z0, z1;

	Re(z0) = 1.0;
	Im(z0) = 0.5;
	z1 = z0;
	
	CosZ(z0);
	SinZ(z1);
	printc("", z0, "\n");
	printc("", z1, "\n");
	
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
