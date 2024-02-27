#include <stdio.h>
#include <math.h>

#include "complex.h"

double _M_PI;
double _M_2PI, _M_LOG10;

int
main()
{
	complex_t a, phasor, cacc, erracc, t0;
	double ph, phacc;
	int j;
	
	_M_PI = 4.0 * atan(1.0);
	_M_2PI = 8.0 * atan(1.0);
	_M_LOG10 = log(10.0);

	ph = _M_2PI / 960.0;
	expiphi(a, _M_PI / 4.0);

	for (j = 0, phacc = 0.0,
			Re(cacc) = 0.0, Im(cacc) = 0.0,
			Re(erracc) = 0.0, Im(erracc) = 0.0; ; j++) {
		expiphi(phasor, phacc);
		printf("%d %.16lf ", j, Arg(phasor));
		madd(cacc, a, phasor);
		printc(cacc);
		printf(" ");
		LogZ2(t0, cacc);
		expZ(t0);
		printc(t0);
		printf(" ");
		sub2(t0, cacc);
		add2(erracc, t0);
#ifdef ERRACC
		printc(erracc);
		printf("\n");
#else
		printf("%.16lf\n", log(mag(erracc)) / _M_LOG10);
#endif
		phacc += ph;
	}

	return 0;
}