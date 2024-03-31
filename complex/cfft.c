#include <math.h>

#include "complex.h"
#include "cfft.h"

static void
rearrange(N, data)
	long N;
	complex_t *data;
{
	long position, mask;
	long target = 0l;
	complex_t t0;
	
	for (position = 0l; position < N; position++) {
		if (target > position) {
			t0 = data[target];
			data[target] = data[position];
			data[position] = t0;
		}
		mask = N;
		while (target & (mask >>= 1))
			target &= ~mask;
		target |= mask;
	}
	
	return;
}

static void
rearrange2(N, outdata, indata)
	long N;
	complex_t *outdata;
	complex_t *indata;
{
	long position, mask;
	long target = 0l;
	
	for (position = 0l; position < N; position++) {
		outdata[target] = indata[position];
		mask = N;
		while (target & (mask >>= 1))
			target &= ~mask;
		target |= mask;
	}
	
	return;
}

static void
do_fft(N, data, inverse)
	long N;
	complex_t *data;
	int inverse;
{
	long step, jump, group, pair, match;
	complex_t multiplier, factor, product;
	double pi, delta, sine;
	
	pi = inverse ? M_PI : -M_PI;
	
	for (step = 1l; step < N; step <<= 1) {
		jump = step << 1;
		delta = pi / (double)step;
		sine = sin(delta * 0.5);
		mkC(multiplier, -2.0 * sine * sine, sin(delta));
		mk1(factor);
		
		for (group = 0u; group < step; group++) {
			for (pair = group; pair < N; pair += jump) {
				match = pair + step;
				mul3(product, factor, data[match]);
				sub3(data[match], data[pair], product);
				add2(data[pair], product);
			}
			madd(factor, factor, multiplier);
		}
	}
	
	return;
}

static void
scaleout(N, data)
	long N;
	complex_t *data;
{
	long position;
	double factor;
	
	factor = 1.0 / (double)N;
	
	for (position = 0l; position < N; position++)
		scale2(data[position], factor);
	
	return;
}

int
fft2(N, outdata, indata)
	long N;
	complex_t *outdata;
	complex_t *indata;
{
	if ((N < 1l) || (N & (N - 1l))) return 0;
	
	rearrange2(N, outdata, indata);
	do_fft(N, outdata, 0);
	
	return 1;
}

int
fft(N, data)
	long N;
	complex_t *data;
{
	if ((N < 1l) || (N & (N - 1l))) return 0;
	
	rearrange(N, data);
	do_fft(N, data, 0);
	
	return 1;
}
	
int
ifft2(N, outdata, indata, scale)
	long N;
	complex_t *outdata;
	complex_t *indata;
	int scale;
{
	if ((N < 1l) || (N & (N - 1l))) return 0;
	
	rearrange2(N, outdata, indata);
	do_fft(N, outdata, 1);
	if (scale) scaleout(N, outdata);
	
	return 1;
}
	
int
ifft(N, data, scale)
	long N;
	complex_t *data;
	int scale;
{
	if ((N < 1l) || (N & (N - 1l))) return 0;
	
	rearrange(N, data);
	do_fft(N, data, 1);
	if (scale) scaleout(N, data);
	
	return 1;
}

