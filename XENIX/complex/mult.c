#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"
#include "cmatrix.h"

#define SIZE 400l
#define SIZESQ (SIZE * SIZE)

complex_t *c0, *c1, *a, *b, *bt;
complex_t *marr;

void compare(long, complex_t *, complex_t *);

int
main()
{
	long i, j;
	int cmp;
	
	printf("init\n");
	
	marr = (complex_t *)calloc(5l * SIZESQ, sizeof (complex_t));
	if (!marr) {
		fprintf(stderr, "cannot alloc marr\n");
		goto bad0;
	} else printf("marr %16lx % 10ld\n",
		(unsigned long)marr, 5l * SIZESQ * sizeof (complex_t));

	a = marr + 0l * SIZESQ;
	c0 = marr + 1l * SIZESQ;
	b = marr + 2l * SIZESQ;
	c1 = marr + 3l * SIZESQ;
	bt = marr + 4l * SIZESQ;
	printf("a %16lx | b %16lx | bt %16lx | c0 %16lx | c1 %16lx\n",
		(unsigned long)a, (unsigned long)b, (unsigned long)bt,
		(unsigned long)c0, (unsigned long)c1);

	for (i = 0l; i < SIZE; i++)
		for (j = 0l; j < SIZE; j++)
			mkZ(MIJ(SIZE, a, i, j),
				(double)((3l * i) % SIZE) / 97.0,
				(double)((7l * j) % SIZE) / 91.0);

	contransAB(SIZE, b, a);

	printf("begin\n");

	printf("c0 %16lx\n", (unsigned long)c0);
	mulCAB(SIZE, c0, a, b);
	
	printf("c1 %16lx\n", (unsigned long)c1);
#ifdef TRANS
	transAB(SIZE, bt, b);
#else
	copyAB(SIZE, bt, b);
	transA(SIZE, bt);
#endif
	mulCABtrans(SIZE, c1, a, bt);

	cmp = memcmp((void *)c0, (void *)c1, SIZESQ * sizeof (complex_t));
	printf("compare %d\n", cmp);

	printf("end\n");

	if (cmp) compare(SIZESQ, c0, c1);
	
	free((void *)marr);
	
bad0:
	return 0;
}

void
compare(size, c0, c1)
	long size;
	complex_t *c0;
	complex_t *c1;
{
	long i, j;
	int cmp;
	unsigned char *pb0, b0;
	unsigned char *pb1, b1;
	
	for (i = 0l; i < size; i++)
		if ((cmp = memcmp(
				(void *)(c0 + i),
				(void *)(c1 + i),
				sizeof (complex_t)))) {
			printf("%ld %d %16lx %16lx\n", i, cmp,
				(unsigned long)c0 + i, (unsigned long)c1 + i);
			printZ("", c0[i], " | ");
			printZ("", c1[i], "\n");
			for (j = 0l; j < sizeof (complex_t); j++) {
				pb0 = (unsigned char *)(c0 + i) + j;
				b0 = *pb0;
				pb1 = (unsigned char *)(c1 + i) + j;
				b1 = *pb1;
				printf("\t%16lx %02x | %16lx %02x%s",
					(unsigned long)pb0, (unsigned int)b0,
					(unsigned long)pb1, (unsigned int)b1,
					(b0 == b1) ? "\n" : " ***\n");
			}
		}
	
	return;
}

