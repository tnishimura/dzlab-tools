#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include "countvec.h"
 
int main(int argc, char * argv[]) {
    countvec_t * v = countvec_new(20, 1, 1);
    countvec_dump(v);
    countvec_increment_range(v, 2, 4, 3.3);
    countvec_increment_range(v, 200, 4, 3.3);
    countvec_multiply_range(v, 3, 5, 3.3);
    countvec_dump(v);
    countvec_destroy(&v);
 
    return 0;
}
