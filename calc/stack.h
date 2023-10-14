#ifndef _STACK_
#define _STACK_

#include <stdlib.h>
#include <stdint.h>

#include "datum.h"

typedef struct stack stack_t;
typedef int (*stack_iter_t)(void *, data_t);

stack_t *stack_constructor(size_t size);
void stack_destructor(stack_t *s);
int stack_push(stack_t *s, data_t val);
int stack_pop(stack_t *s, data_t *val);
int stack_peek(stack_t *s, data_t *val);
int stack_iter_peek(stack_t *s, stack_iter_t iter, void *_ctx);

#endif

