#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "sm.h"
#include "stack.h"

#define STATE_START 0
#define STATE_OPERATOR 1
#define STATE_NUMBER 2
#define STATE_FINAL 3

#define EVENT_DELIM 0
#define EVENT_DIGIT 1
#define EVENT_OP 2
#define EVENT_TERM 3

#define STACK_SIZE 1024ul

struct ctx {
	int base;
	int acc;
	int acc_valid;
	char digit[2];
	stack_t *s;
};

static char *hexp = "PT %0.8x\n";
static char *decp = "PT %d\n";
static char *hexP = "PA %0.8x\n";
static char *decP = "PA %d\n";
static char *hexv = "VT %0.8x\n";
static char *decv = "VT %d\n";
static char *hexV = "VA %0.8x\n";
static char *decV = "VA %d\n";

static int accumulate(void *_ctx, delta_t *delta) {
	struct ctx *ctx = (struct ctx *)_ctx;
	
	if (!ctx->acc_valid) {
		ctx->acc_valid = 1;
		ctx->acc = 0;
	}
	
	ctx->acc = ctx->acc * ctx->base + (int)strtol(ctx->digit, NULL, ctx->base);
	
	return 1;
}

static int dump_pop(void *_ctx, delta_t *delta) {
	struct ctx *ctx = (struct ctx *)_ctx;
	data_t d;
	int r;
	
	r = stack_pop(ctx->s, &d);
	if (!r) {
		if (delta->event != EVENT_TERM) printf("stack underflow\n");
	} else while (r > 0) {
		if (d.type == INT) printf((ctx->base == 16) ? hexV : decV, d.datum.i);
		r = stack_pop(ctx->s, &d);
	}
	
	return 1;
}

static int dump_peek(void *_ctx, data_t d) {
	struct ctx *ctx = (struct ctx *)_ctx;
	
	if (d.type == INT) printf((ctx->base == 16) ? hexP : decP, d.datum.i);
	
	return 1;
}

static int operator(void *_ctx, delta_t *delta) {
	struct ctx *ctx = (struct ctx *)_ctx;
	data_t d0, d1;
	
	switch (ctx->digit[0]) {
	case 'h':
		ctx->base = 10;
		break;
	case 'H':
		ctx->base = 16;
		break;
	case 'p':
		if (!stack_peek(ctx->s, &d0)) printf("stack underflow\n");
		else if (d0.type == INT) printf((ctx->base == 16) ? hexp : decp, d0.datum.i);
		break;
	case 'P':
		if (!stack_iter_peek(ctx->s, dump_peek, ctx)) printf("stack underflow\n");
		break;
	case '.':
	case 'v':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (d0.type == INT) printf((ctx->base == 16) ? hexv : decv, d0.datum.i);
		break;
	case 'V':
		(void)dump_pop(_ctx, delta);
		break;
	case 'u':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else {
			(void)stack_push(ctx->s, d0);
			if (!stack_push(ctx->s, d0)) printf("stack overflow\n");
		}
		break;
	case '+':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i += d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	case '-':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i -= d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	case '*':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i *= d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	case '/':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i /= d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	case '%':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i %= d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	case '&':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i &= d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	case '|':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i |= d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	case '~':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (d0.type == INT) {
			d0.datum.i = ~d0.datum.i;
			(void)stack_push(ctx->s, d0);
		} else return UNDEF;
		break;
	case '^':
		if (!stack_pop(ctx->s, &d0)) printf("stack underflow\n");
		else if (!stack_pop(ctx->s, &d1)) {
			(void)stack_push(ctx->s, d0);
			printf("stack underflow\n");
		} else {
			if ((!(d0.type == INT)) || (!(d1.type == INT)))
				return UNDEF;
			d0.datum.i ^= d1.datum.i;
			(void)stack_push(ctx->s, d0);
		}
		break;
	default:
		return UNDEF;
	}
	
	return 1;
}

static int push_acc(void *_ctx, delta_t *delta) {
	struct ctx *ctx = (struct ctx *)_ctx;
	data_t d;
	
	ctx->acc_valid = 0;
	d.type = INT;
	d.datum.i = ctx->acc;
	
	if (!stack_push(ctx->s, d)) printf("stack overflow\n");
	
	if (delta->event == EVENT_OP) return operator(_ctx, delta);
	else return 1;
}

static delta_t deltas[] = {
	{ STATE_START, EVENT_DELIM, STATE_START, NULL, NULL},
	{ STATE_START, EVENT_DIGIT, STATE_NUMBER, NULL, accumulate },
	{ STATE_START, EVENT_OP, STATE_OPERATOR, NULL, operator },
	
	{ STATE_OPERATOR, EVENT_DELIM, STATE_START, NULL, NULL },
	{ STATE_OPERATOR, EVENT_OP, STATE_OPERATOR, NULL, operator },
	{ STATE_OPERATOR, EVENT_DIGIT, STATE_NUMBER, NULL, accumulate },
	
	{ STATE_NUMBER, EVENT_DELIM, STATE_START, NULL, push_acc },
	{ STATE_NUMBER, EVENT_DIGIT, STATE_NUMBER, NULL, accumulate },
	{ STATE_NUMBER, EVENT_OP, STATE_OPERATOR, NULL, push_acc },
	
	{ ANY, EVENT_TERM, STATE_FINAL, NULL, dump_pop },
	
	{ UNDEF, UNDEF, UNDEF, NULL, NULL }
};

int main(void) {
	FILE *stream = stdin;
	struct ctx c = { 10, 0, 0, { '\0', '\0' }, NULL };
	int input;
	state_t s;
	int ret;
	
	if (!(c.s = stack_constructor(STACK_SIZE))) return -1;
	
	state_init(&s, STATE_START, STATE_FINAL, UNDEF, deltas, &c);
	
	while (1) {
		input = fgetc(stream);
		if ((input == EOF) || ((char)input == 'q')) {
			if (state_exec(&s, EVENT_TERM) <= 0) break;
		} else if (isxdigit(input)) {
			c.digit[0] = (char)input;
			if (state_exec(&s, EVENT_DIGIT) <= 0) break;
		} else if (
				((char)input == 'h') || ((char)input == 'H') ||
				((char)input == 'p') || ((char)input == 'P') ||
				((char)input == '.') || ((char)input == 'v') || ((char)input == 'V') ||
				((char)input == 'u') ||
				((char)input == '+') ||
				((char)input == '-') ||
				((char)input == '*') ||
				((char)input == '/') ||
				((char)input == '%') ||
				((char)input == '&') ||
				((char)input == '|') ||
				((char)input == '~') ||
				((char)input == '^')
		) {
			c.digit[0] = (char)input;
			if (state_exec(&s, EVENT_OP) <= 0) break;
		} else {
			if (state_exec(&s, EVENT_DELIM) <= 0) break;
		}
	}
	
	stack_destructor(c.s);
	
	return 0;
}

