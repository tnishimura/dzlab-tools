#ifndef _CDZUTILS_H
#define _CDZUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct{
    int base;
    int size;
    int first;
    int last;
    double * data;
} countvec_t;

// {con,de}structor
countvec_t * countvec_new(int size, int base, double initial);
void countvec_destroy(countvec_t ** v);

// other
void countvec_dump(countvec_t * v);
int countvec_is_in_range(countvec_t * v, int from, int to);
int sizeof_int();
int sizeof_double();

// getters
int countvec_pack_int(countvec_t * v, void * vbuffer, int from, int to);
int countvec_pack_double(countvec_t * v, void * vbuffer, int from, int to);
double countvec_get(countvec_t * v, int pos);
int countvec_first(countvec_t * v);
int countvec_last (countvec_t * v);

// mutators
int countvec_increment_range(countvec_t * v, int from, int to, double value);
int countvec_multiply_range(countvec_t * v, int from, int to, double value);
int countvec_set_range(countvec_t * v, int from, int to, double value);

#endif
