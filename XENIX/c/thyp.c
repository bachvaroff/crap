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
	
	CoshZ(z0);
	SinhZ(z1);
	printZ("", z0, "\n");
	printZ("", z1, "\n");
	mulZ2(z0, z0);
	mulZ2(z1, z1);
	printZ("", z0, "\n");
	printZ("", z1, "\n");
	subZ3(u, z0, z1);
	printZ("", u, "\n");
	
	LogZ(z0);
	printZ("", z0, "\n");
	
	expZ(z0);
	printZ("", z0, "\n");
	
	LogZ(z1);
	printZ("", z1, "\n");
	
	expZ(z1);
	printZ("", z1, "\n");
		
	mkZ(z0, log(1.0 + sqrt(5.0)) - M_LN2, 0.0);
	SinhZ2(z1, z0);
	printZ("", z1, "\n");
	CoshZ2(z1, z0);
	printZ("", z1, "\n");
	TanhZ2(z1, z0);
	printZ("", z1, "\n");
	TanhZ(z0);
	printZ("", z0, "\n");

	mkZ0(z0);
	TanhZ(z0);
	printZ("", z0, "\n");
	
	mkZ1(z0);
	TanhZ(z0);
	printZ("", z0, "\n");

	mkZi(z0);
	TanhZ(z0);
	printZ("", z0, "\n");

	mkZ(z0, 1.0, 1.0);
	TanhZ(z0);
	printZ("", z0, "\n");
	
	return 0;
}

