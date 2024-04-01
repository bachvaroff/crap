#include <math.h>

#include "complex.h"
#include "cfft.h"

static void
rearrange(N, outdata, indata)
	long N;
	complex_t *outdata;
	complex_t *indata;
{
	long position, target, mask;
	
	for (position = 0l, target = 0l; position < N; position++) {
		if (outdata == indata) {
			if (target > position)
				xchgZ(outdata[target], outdata[position]);
		} else outdata[target] = indata[position];
		for (mask = N >> 1; target & mask; mask >>= 1)
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
	double pi, theta;
	
	pi = inverse ? M_PI : -M_PI;
	
	for (step = 1l; step < N; step <<= 1) {
		jump = step << 1;
		theta = pi / (double)step;
		mkZ1(factor);
		expiphiZ(multiplier, theta);
		subZ2(multiplier, factor);
		
		for (group = 0u; group < step; group++) {
			for (pair = group; pair < N; pair += jump) {
				match = pair + step;
				mulZ3(product, factor, data[match]);
				subZ3(data[match], data[pair], product);
				addZ2(data[pair], product);
			}
			maddZ(factor, factor, multiplier);
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
		scaleZ2(data[position], factor);
	
	return;
}

int
fft2(N, outdata, indata)
	long N;
	complex_t *outdata;
	complex_t *indata;
{
	if ((N < 1l) || (N & (N - 1l))) return 0;
	
	rearrange(N, outdata, indata);
	do_fft(N, outdata, 0);
	
	return 1;
}

int
fft(N, data)
	long N;
	complex_t *data;
{
	if ((N < 1l) || (N & (N - 1l))) return 0;
	
	rearrange(N, data, data);
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
	
	rearrange(N, outdata, indata);
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
	
	rearrange(N, data, data);
	do_fft(N, data, 1);
	if (scale) scaleout(N, data);
	
	return 1;
}

