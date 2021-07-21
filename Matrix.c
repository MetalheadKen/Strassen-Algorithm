#include <assert.h>
#include <stdlib.h>
#include "Matrix.h"

static inline int **matrix_allocate_value(uint32_t row, uint32_t column)
{
    int **matrix = (int **) calloc(row + row * column, sizeof(**matrix));
    int *matrixTemp = (int *) (matrix + row);

    for (uint32_t i = 0; i < row; i++) {
        matrix[i] = matrixTemp;
        matrixTemp += column;
    }

    return matrix;
}

static inline void Matrix_Free(Matrix **ptr)
{
    assert(*ptr && (*ptr)->values);
    free((*ptr)->values);
    free(*ptr);
    *ptr = 0 ;
}

#define matrix_autofree __attribute__((cleanup(Matrix_Free)))

static inline Matrix * Matrix_Addition(Matrix * res, const Matrix * a, const Matrix * b)
{
    assert( res && a && b);
    for (uint32_t i = 0; i < res->row; i++)
        for (uint32_t j = 0; j < res->column; j++)
            res->values[i][j] = a->values[i][j] + b->values[i][j];

    return res;
}

static inline Matrix * Matrix_Subtract(Matrix * res, const Matrix * a, const Matrix * b)
{
    assert( res && a && b);
    for (uint32_t i = 0; i < res->row; i++)
        for (uint32_t j = 0; j < res->column; j++)
            res->values[i][j] = a->values[i][j] - b->values[i][j];

    return res;
}

static inline Matrix * Matrix_Multiply(Matrix * res, const Matrix * a, const Matrix * b)
{
    assert( res && a && b);
    int m1 = (a->values[0][0] + a->values[1][1]) * (b->values[0][0] + b->values[1][1]);
    int m2 = (a->values[1][0] + a->values[1][1]) *  b->values[0][0];
    int m3 = (b->values[0][1] - b->values[1][1]) *  a->values[0][0];
    int m4 = (b->values[1][0] - b->values[0][0]) *  a->values[1][1];
    int m5 = (a->values[0][0] + a->values[0][1]) *  b->values[1][1];
    int m6 = (a->values[1][0] - a->values[0][0]) * (b->values[0][0] + b->values[0][1]);
    int m7 = (a->values[0][1] - a->values[1][1]) * (b->values[1][0] + b->values[1][1]);

    res->values[0][0] = m1 + m4 - m5 + m7;
    res->values[0][1] = m3 + m5;
    res->values[1][0] = m2 + m4;
    res->values[1][1] = m1 - m2 + m3 + m6;
  
    return res;
}

static inline Matrix * Matrix_Initializer( uint32_t row, uint32_t column)
{
    Matrix *matrix = calloc(1, sizeof(Matrix));
    matrix->row    = row;
    matrix->column = column;
    matrix->values = matrix_allocate_value(row, column);
    return matrix ;
}
