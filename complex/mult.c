#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"
#include "cmatrix.h"

#define SIZE 400u
#define SIZESQ (SIZE * SIZE)

complex_t *c0, *c1, *a, *b, *bt;
complex_t *marr;

void compare(size_t, complex_t *, complex_t *);

int
main()
{
	size_t i, j;
	int cmp;
	
	printf("init\n");
	
	marr = (complex_t *)calloc(5u * SIZESQ, sizeof (complex_t));
	if (!marr) {
		fprintf(stderr, "cannot alloc marr\n");
		goto bad0;
	} else printf("marr %08x % 10u\n", marr, 5u * SIZESQ * sizeof (complex_t));

	a = marr + 0u * SIZESQ;
	c0 = marr + 1u * SIZESQ;
	b = marr + 2u * SIZESQ;
	c1 = marr + 3u * SIZESQ;
	bt = marr + 4u * SIZESQ;
	printf("a %08x b %08x bt %08x c0 %08x c1 %08x\n", a, b, c0, c1);

	for (i = 0u; i < SIZE; i++)
		for (j = 0u; j < SIZE; j++) {
			Re(IDX(a, SIZE, i, j)) = (double)((3u * i) % SIZE) / 97.0;
			Im(IDX(a, SIZE, i, j)) = (double)((7u * j) % SIZE) / 91.0;
		}

	contransAB(SIZE, b, a);

	printf("begin\n");

	printf("c0 %08x\n", c0);
	mulCAB(SIZE, c0, a, b);
	
	printf("c1 %08x\n", c1);
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
	size_t size;
	complex_t *c0;
	complex_t *c1;
{
	size_t i, j;
	int cmp;
	unsigned char *pb0, b0;
	unsigned char *pb1, b1;
	
	for (i = 0u; i < size; i++)
		if ((cmp = memcmp((void *)(c0 + i), (void *)(c1 + i), sizeof (complex_t)))) {
			printf("%u %d %08x %08x\n", i, cmp, c0 + i, c1 + i);
			printc("", c0[i], " | ");
			printc("", c1[i], "\n");
			for (j = 0u; j < sizeof (complex_t); j++) {
				pb0 = (unsigned char *)(c0 + i) + j;
				b0 = *pb0;
				pb1 = (unsigned char *)(c1 + i) + j;
				b1 = *pb1;
				printf("\t%08x %02x | %08x %02x%s",
						pb0, (unsigned int)b0,
						pb1, (unsigned int)b1,
						(b0 == b1) ? "\n" : " ***\n");
			}
		}
	
	return;
}
