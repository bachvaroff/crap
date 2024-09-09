#include <stdio.h>
#include <math.h>

#include "complex.h"

#ifndef STEPS
#define STEPS 4096
#endif

int
main()
{
	complex_t a, phasor, cacc, erracc, t0;
	double ph, phacc, t1;
	int j, jbr, jph;
	
	ph = 2.0 * M_PI / (double)STEPS;
	
	mkZi(a);
	mkZ0(cacc);
	mkZ0(erracc);
	for (j = 0, phacc = 0.0; ; j++) {
		jbr = j / STEPS;
		jph = j % STEPS;
		expiphiZ(phasor, phacc);
		printf("%d %d %.16lf %.16lf ", jbr, jph,
				(double)(jph) * ph, ArgZ(phasor));
		maddZ(cacc, a, phasor);
		printZ("", cacc, " ");
		logZ2(t0, cacc, jbr);
		printZ("", t0, " ");
		expZ(t0);
		printZ("", t0, " ");
		subZ2(t0, cacc);
		addZ2(erracc, t0);
		t1 = magZ(erracc);
		printf("%.16lf\n", log(t1) / M_LN10);
		phacc += ph;
	}

	return 0;
}
