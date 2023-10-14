#ifndef _DATUM_
#define _DATUM_

typedef enum datum_type {
	CHAR,
	SHORT,
	INT,
	LONG,
	POINTER
} datum_type_t;

typedef union datum {
	char c;
	short s;
	int i;
	long l;
	void *p;
} datum_t;

typedef struct data {
	datum_type_t type;
	datum_t datum;
} data_t;

#endif

