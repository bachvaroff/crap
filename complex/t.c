#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

int
main()
{
	complex_t z0, z1, u;

	Re(z0) = 3.0;
	Im(z0) = 0.0;
	z1 = z0;
	
	CoshZ(z0);
	SinhZ(z1);
	printc("", z0, "\n");
	printc("", z1, "\n");
	mul2(z0, z0);
	mul2(z1, z1);
	printc("", z0, "\n");
	printc("", z1, "\n");
	sub3(u, z0, z1);
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
