#ifndef _CFFT_H
#define _CFFT_H

#include "complex.h"

void scaledata(long, complex_t *);
int fft(long, complex_t *);
int fft2(long, complex_t *, complex_t *);
int dft2(long, complex_t *, complex_t *);
int ifft(long, complex_t *, int);
int ifft2(long, complex_t *, complex_t *, int);
int idft2(long, complex_t *, complex_t *, int);

#endif

