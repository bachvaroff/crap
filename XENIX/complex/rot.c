#include <stdio.h>
#include <math.h>

#include "complex.h"

#ifndef STEPS
#define STEPS 4096.0
#endif

int
main()
{
	complex_t a, phasor, cacc, erracc, t0;
	double ph, phacc;
	int j;
	
	ph = 2.0 * M_PI / STEPS;
	expiphiZ(a, M_PI_4);
	
	mkZ0(cacc);
	mkZ0(erracc);
	for (j = 0, phacc = 0.0; ; j++) {
		expiphiZ(phasor, phacc);
		printf("%d %.16lf ", j, ArgZ(phasor));
		maddZ(cacc, a, phasor);
		printZ("", cacc, " ");
		LogZ2(t0, cacc);
		expZ(t0);
		printZ("", t0, " ");
		subZ2(t0, cacc);
		printZ("", t0, " ");
		addZ2(erracc, t0);
#ifdef ERRACC
		printZ("", erracc, "\n");
#else
		printf("%.16lf\n", log(magZ(erracc)) / M_LN10);
#endif
		phacc += ph;
	}

	return 0;
}
