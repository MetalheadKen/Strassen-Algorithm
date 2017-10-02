#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdint.h>

#define MATRIX_INITIALIZER(X, ROW, COLUMN) \
        Matrix_Initializer(&(X) ,(ROW), (COLUMN));

enum { NAIVE_ARITHMETIC, };

typedef struct _Matrix {
    uint32_t row;
    uint32_t column;
    int **values;
    
    int **(*New)(uint32_t row, uint32_t column);
    void  (*Delete)(void *);
} Matrix;

typedef struct _Matrix_Arith {
    Matrix (*Addition)(Matrix, const Matrix, const Matrix);
    Matrix (*Subtract)(Matrix, const Matrix, const Matrix);
    Matrix (*Multiply)(Matrix, const Matrix, const Matrix);
} Matrix_Arith;

void Matrix_Initializer(Matrix *, uint32_t, uint32_t);

extern Matrix_Arith Naive_Matrix_Arith;

#endif /* MATRIX_H_ */
