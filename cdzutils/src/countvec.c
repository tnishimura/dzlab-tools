#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "countvec.h"

// constructor
countvec_t * countvec_new(int size, int base, double initial){
    countvec_t * v = malloc(sizeof(countvec_t));
    assert(v);

    v->data = malloc(sizeof(double) * (size + 1));
    assert(v->data);
    v->size  = size;
    v->base  = base;
    v->first = base;
    v->last  = v->size - 1 + v->base;
    for (int i = v->first; i <= v->last; i++) {
        v->data[i] = initial;
    }
    return v;
}

// deconstructor
void countvec_destroy(countvec_t ** v){
    printf("countvec_destroy called\n");
    assert(v); assert(*v);
    free((*v)->data);
    free(*v);
    *v = NULL;
    printf("countvec_destroy finished\n");
}

void countvec_dump(countvec_t * v){
    assert(v);

    for (int i = v->first; i <= v->last; i++) {
        printf("%d: %f\n", i, v->data[i]);
    }
}

int countvec_is_in_range(countvec_t * v, int from, int to){
    assert(v);
    return v->first <= from && from <= to && to <= v->last;
}

// assume buffer is of size (to - from + 1) * sizeof(double)
// faster to use this and unpack() in perl, than calling countvec_get multiple times
int countvec_pack_double(countvec_t * v, void * vbuffer, int from, int to){
    assert(v);
    if (countvec_is_in_range(v, from, to)){
        double * dbuffer = vbuffer;
        for (int i = from; i <= to; i++) {
            dbuffer[i - from] = v->data[i];
        }
        return 1;
    }
    return 0;
}

// same as above but round to int
int countvec_pack_int(countvec_t * v, void * vbuffer, int from, int to){
    assert(v);
    if (countvec_is_in_range(v, from, to)){
        int * ibuffer = vbuffer;
        for (int i = from; i <= to; i++) {
            ibuffer[i - from] = v->data[i] >= 0 ? (int)(v->data[i] + 0.5) : (int)(v->data[i] - 0.5);
        }
        return 1;
    }
    return 0;
}

double countvec_get(countvec_t * v, int pos){
    assert(v);
    if (countvec_is_in_range(v, pos, pos)){
        return v->data[pos];
    }
    return NAN;
}

int countvec_first(countvec_t * v){ assert(v); return v->first; }
int countvec_last (countvec_t * v){ assert(v); return v->last; }

///////////////////////////////////////////////////////////
// mutators

int countvec_increment_range(countvec_t * v, int from, int to, double value){
    assert(v);
    if (! countvec_is_in_range(v, from ,to)){
        return 0;
    }
    for (int i = from; i <= to; i++) {
        v->data[i] += value;
    }
    return 1;
}

int countvec_set_range(countvec_t * v, int from, int to, double value){
    assert(v);
    if (! countvec_is_in_range(v, from ,to)){
        return 0;
    }
    for (int i = from; i <= to; i++) {
        v->data[i] = value;
    }
    return 1;
}

int countvec_multiply_range(countvec_t * v, int from, int to, double value){
    assert(v);
    if (! countvec_is_in_range(v, from ,to)){
        return 0;
    }
    for (int i = from; i <= to; i++) {
        v->data[i] *= value;
    }
    return 1;
}

int sizeof_int(){
    return sizeof(int);
}

int sizeof_double(){
    return sizeof(double);
}

