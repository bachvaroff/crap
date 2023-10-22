#include <stdlib.h>
#include "stack.h"

typedef struct stack {
	data_t *data;
	size_t N;
	ssize_t p;
} stack_t;

stack_t *stack_constructor(size_t size) {
	stack_t *s;
	
	s = (stack_t *)malloc(sizeof (stack_t));
	if (!s) goto bad0;
	
	s->data = (data_t *)malloc(size * sizeof (data_t));
	if (!s->data) goto bad1;
	
	s->N = size;
	s->p = -1l;
	
	return s;
	
bad1:
	free(s);
bad0:
	return NULL;
}

void stack_destructor(stack_t *s) {
	if (!s) return;
	if (s->data) free(s->data);
	s->data = NULL;
	free(s);
	return;
}

int stack_push(stack_t *s, data_t val) {
	if (!s->data) return -1;
	if (s->p == (s->N - 1)) return 0;
	s->p++;
	s->data[s->p] = val;
	return 1;
}

int stack_pop(stack_t *s, data_t *val) {
	if (!s->data) return -1;
	if (s->p < 0) return 0;
	*val = s->data[s->p];
	s->p--;
	return 1;
}

int stack_peek(stack_t *s, data_t *val) {
	if (!s->data) return -1;
	if (s->p < 0) return 0;
	*val = s->data[s->p];
	return 1;
}

int stack_iter_peek(stack_t *s, stack_iter_t iter, void *_ctx) {
	int j, r;
	
	if (!s->data) return -1;
	if (s->p < 0) return 0;
	
	for (j = s->p, r = 0; j >= 0; j--)
		if ((r = iter(_ctx, s->data[j])) <= 0) break;
	
	return r;
}

