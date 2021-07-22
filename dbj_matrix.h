#ifndef MATRIX_H_
#define MATRIX_H_

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct _Matrix
{
    unsigned rows;
    unsigned cols;
    int *values;
} Matrix;

#define MXY_NOT_POINTER_(M_, R_, C_) ((int(*)[M_.rows])(M_.values))[R_][C_]

// mx as a pointer
#define MXY(M_, R_, C_) MXY_NOT_POINTER_(((Matrix)*M_), R_, C_)

#define FOR(C_, N_) for (unsigned C_ = 0; C_ < N_; ++C_)

/*
    Matrix mx = { 3,2, calloc(3 * 2, sizeof(int))  } ;

    MXY(mx,1,1) = 42;

    assert( MXY(mx,1,1) == 42) ;
*/

static inline void matrix_free(Matrix **ptr)
{
    assert(*ptr && (*ptr)->values);
    free((*ptr)->values);
    free(*ptr);
    *ptr = 0;
}

#define matrix_autofree __attribute__((cleanup(matrix_free)))

static inline Matrix *matrix_addition(Matrix *res, const Matrix *a, const Matrix *b)
{
    assert(res && a && b);
    FOR(i, res->rows)
    FOR(j, res->cols)
    // res->values[i][j] = a->values[i][j] + b->values[i][j];
    MXY(res, i, j) = MXY(a, i, j) + MXY(b, i, j);

    return res;
}

static inline Matrix *matrix_subtraction(Matrix *res, const Matrix *a, const Matrix *b)
{
    assert(res && a && b);
    FOR(i, res->rows)
    FOR(j, res->cols)
    MXY(res, i, j) = MXY(a, i, j) - MXY(b, i, j);

    return res;
}

static inline int matrix_equal(const Matrix *a, const Matrix *b)
{
    assert(a && b);

    if (a->rows != b->rows)
        return 0;
    if (a->cols != b->cols)
        return 0;

    FOR(i, a->rows)
    FOR(j, a->cols)
    if (MXY(a, i, j) != MXY(b, i, j))
        return 0;

    return 1;
}

static Matrix *matrix_ijk_matmul(Matrix *C, const Matrix *A, const Matrix *B)
{
    assert(A && B && C);
    FOR(i, A->rows)
    {
        FOR(j, B->cols)
        {
            MXY(C, i, j) = 0;
            FOR(k, A->cols)
            {
                MXY(C, i, j) += MXY(A, i, k) * MXY(B, k, j);
            }
        }
    }
    return C;
}

static inline Matrix *matrix_new(uint32_t rows, uint32_t cols)
{
    Matrix *matrix = calloc(1, sizeof(Matrix));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->values = calloc(rows * cols, sizeof(int));
    return matrix;
}

static void matrix_print(const char *prompt, Matrix *mx)
{
    assert(prompt && mx);
    printf("%s\n", prompt);
    FOR(i, mx->rows)
    {
        FOR(j, mx->cols)
        {
            printf(" %5d ", MXY(mx, i, j));
        }
        printf("\n");
    }
}

static void inline matrix_foreach(Matrix *mx, Matrix *(*callback)(Matrix *, unsigned, unsigned))
{
    FOR(i, mx->rows)
    {
        FOR(j, mx->cols)
        {
            callback(mx, i, j);
        }
    }
}

#endif /* MATRIX_H_ */
