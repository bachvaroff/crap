#include <stdio.h>
#include <math.h>

#include "complex.h"
#include "cfft.h"

int
main()
{
	complex_t arr[1024];
	complex_t oarr[1024];
	double w;
	long i, j;
	
/* 1024 samples of a square wave @ 8000sps, frequency 62.5Hz */
/*
	for (j = 0l; j < 16l; j++)
		for (i = 0l; i < 64l; i++)
			mkZ(arr[j * 64l + i], ((j & 1l) ? 1.0 : -1.0) / 32768.0, 0.0);
*/
	
/* 1024 samples of a "sawtooth" @ 8000sps, frequency 125Hz */
	for (j = 0l; j < 16l; j++)
		for (i = 0l; i < 64l; i++)
			mkZ(arr[j * 64l + i], (-32.0 + (double)i) / 32768.0, 0.0);
	
/* DTMF * @ 8000sps */
/*
	for (j = 0l; j < 1024l; j++) {
		w = 2.0 * M_PI * (double)j / 8000.0;
		mkZ(arr[j], sin(941.0 * w) + sin(1209.0 * w), 0.0);
	}
*/
	
	for (j = 0l; j < 1024l; j++) {
		printf("SIG %.16lf ", (double)j / 8000.0);
		printZ("", arr[j], "\n");
	}
	
	fft2(1024l, oarr, arr);
	
	for (j = 0l; j < 1024l; j++) {
		printf("FWDFFT %.16lf ", (double)j * 8000.0 / 1024.0);
		printZ("", oarr[j], " ");
		printf("%.16lf\n", magZ(oarr[j]));
	}
	
	dft2(1024l, oarr, arr);
	
	for (j = 0l; j < 1024l; j++) {
		printf("FWDDFT %.16lf ", (double)j * 8000.0 / 1024.0);
		printZ("", oarr[j], " ");
		printf("%.16lf\n", magZ(oarr[j]));
	}
	
	ifft2(1024l, arr, oarr, 1);
	
	for (j = 0l; j < 1024l; j++) {
		printf("INVFFT %.16lf ", (double)j / 8000.0);
		printZ("", arr[j], "\n");
	}
	
	idft2(1024l, arr, oarr, 1);
	
	for (j = 0l; j < 1024l; j++) {
		printf("INVDFT %.16lf ", (double)j / 8000.0);
		printZ("", arr[j], "\n");
	}
	
	return 0;
}

