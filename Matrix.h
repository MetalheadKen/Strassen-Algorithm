#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdint.h>

enum { NAIVE_ARITHMETIC, };

typedef struct _Matrix {
    uint32_t row;
    uint32_t column;
    int **values;
} Matrix;


#endif /* MATRIX_H_ */
