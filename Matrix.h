#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <stdlib.h>

#define MATRIX_INITIALIZER(X, ROW, COLUMN)          \
    do {                                            \
        (X).row    = ROW;                           \
        (X).column = COLUMN;                        \
        (X).values = NULL;                          \
        (X).New    = Matrix_Allocate;               \
        (X).Delete = free;                          \
                                                    \
        (X).values = (X).New(ROW, COLUMN);          \
    } while (0)

typedef struct _Matrix {
    int row;
    int column;
    int **values;
    
    int **(*New)(int row, int column);
    void  (*Delete)(void *);
} Matrix;

int **Matrix_Allocate(int, int);
void  Matrix_Free(void *);

Matrix Matrix_Addition(Matrix, const Matrix, const Matrix);
Matrix Matrix_Subtract(Matrix, const Matrix, const Matrix);
Matrix Matrix_Multiply(Matrix, const Matrix, const Matrix);

#endif /* MATRIX_H_ */
