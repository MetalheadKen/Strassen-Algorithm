
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dbj_matrix.h"

#define SQUARE_MATRIX_SIDE 64

static inline Matrix *strassen(Matrix *dest, const Matrix *srcA, const Matrix *srcB, int length)
{

    if (length == 2)
        return matrix_ijk_matmul(dest, srcA, srcB);

    int len = length / 2;

    matrix_autofree Matrix *a11 = matrix_new(len, len), *a12 = matrix_new(len, len), *a21 = matrix_new(len, len), *a22 = matrix_new(len, len),
                           *b11 = matrix_new(len, len), *b12 = matrix_new(len, len), *b21 = matrix_new(len, len), *b22 = matrix_new(len, len), *c11 = matrix_new(len, len), *c12 = matrix_new(len, len), *c21 = matrix_new(len, len), *c22 = matrix_new(len, len), *m1 = matrix_new(len, len), *m2 = matrix_new(len, len), *m3 = matrix_new(len, len), *m4 = matrix_new(len, len), *m5 = matrix_new(len, len), *m6 = matrix_new(len, len), *m7 = matrix_new(len, len), *temp1 = matrix_new(len, len), *temp2 = matrix_new(len, len);

    /* Divide matrix into four parts */
    FOR(i, len)
    {
        FOR(j, len)
        {
            MXY(a11, i, j) = MXY(srcA, i, j);
            MXY(a12, i, j) = MXY(srcA, i, j + len);
            MXY(a21, i, j) = MXY(srcA, i + len, j);
            MXY(a22, i, j) = MXY(srcA, i + len, j + len);

            MXY(b11, i, j) = MXY(srcB, i, j);
            MXY(b12, i, j) = MXY(srcB, i, j + len);
            MXY(b21, i, j) = MXY(srcB, i + len, j);
            MXY(b22, i, j) = MXY(srcB, i + len, j + len);
        }
    }

    /* Calculate seven formulas of strassen Algorithm */
    strassen(m1, matrix_addition(temp1, a11, a22), matrix_addition(temp2, b11, b22), len);
    strassen(m2, matrix_addition(temp1, a21, a22), b11, len);
    strassen(m3, a11, matrix_subtraction(temp1, b12, b22), len);
    strassen(m4, a22, matrix_subtraction(temp1, b21, b11), len);
    strassen(m5, matrix_addition(temp1, a11, a12), b22, len);
    strassen(m6, matrix_subtraction(temp1, a21, a11), matrix_addition(temp2, b11, b12), len);
    strassen(m7, matrix_subtraction(temp1, a12, a22), matrix_addition(temp2, b21, b22), len);

    /* Merge the answer of matrix dest */
    /* c11 = m1 + m4 - m5 + m7 = m1 + m4 - (m5 - m7) */
    matrix_subtraction(c11, matrix_addition(temp1, m1, m4), matrix_subtraction(temp2, m5, m7));
    matrix_addition(c12, m3, m5);
    matrix_addition(c21, m2, m4);
    matrix_addition(c22, matrix_subtraction(temp1, m1, m2), matrix_addition(temp2, m3, m6));

    /* Store the answer of matrix multiplication */
    FOR(i, len)
    {
        FOR(j, len)
        {
            MXY(dest, i, j) = MXY(c11, i, j);
            MXY(dest, i, j + len) = MXY(c12, i, j);
            MXY(dest, i + len, j) = MXY(c21, i, j);
            MXY(dest, i + len, j + len) = MXY(c22, i, j);
        }
    }

    return dest;
}

/*
--------------------------------------------------------------------------------------------------------------------
callbacks for matrix_for_each
*/
static inline Matrix *ordinal_as_val(Matrix *mx, unsigned i, unsigned j)
{
    MXY(mx, i, j) = (i * mx->cols + j);
    return mx;
}

static inline Matrix *zoro(Matrix *mx, unsigned i, unsigned j)
{
    assert(mx && mx->values);
    MXY(mx, i, j) = 0;
    return mx;
}

/*
--------------------------------------------------------------------------------------------------------------------
*/
/* Check if value is the power of two or not */
#define ISPOW2(V_) (ceil(log2(V_)) == floor(log2(V_)))

int main(int argc, char *argv[])
{
    const unsigned matrix_side_size = SQUARE_MATRIX_SIDE;

    /* Check if dimensions of matrix is the power of two or not */
    if (!ISPOW2(matrix_side_size))
    {
        printf("\n%s\tERROR: square matrix side must be a power of 2. And current size:%d is not.", argv[0], matrix_side_size);
        return 0;
    }

    matrix_autofree Matrix *matrixA = matrix_new(matrix_side_size, matrix_side_size);
    matrix_autofree Matrix *matrixB = matrix_new(matrix_side_size, matrix_side_size);
    matrix_autofree Matrix *matrixC = matrix_new(matrix_side_size, matrix_side_size);
    matrix_autofree Matrix *matrixR = matrix_new(matrix_side_size, matrix_side_size);

    matrix_foreach(matrixA, ordinal_as_val);
    matrix_foreach(matrixB, ordinal_as_val);
    matrix_foreach(matrixC, zoro);
    matrix_foreach(matrixR, zoro);

    printf("\n\nMatrix multiplication: C = A * B\nMultiplying int[%d][%d] * int[%d][%d]\n R is the result of ijk_matmul(A,B)\n\n", SQUARE_MATRIX_SIDE, SQUARE_MATRIX_SIDE, SQUARE_MATRIX_SIDE, SQUARE_MATRIX_SIDE);

    /* Matrix multiplication */
    matrixC = strassen(matrixC, matrixA, matrixB, matrix_side_size);

    matrixR = matrix_ijk_matmul(matrixR, matrixA, matrixB);

    if (SQUARE_MATRIX_SIDE < 9)
    {
        matrix_print("Matrix A:", matrixA);
        matrix_print("Matrix B:", matrixB);
        matrix_print("Matrix C:", matrixC);
        matrix_print("Matrix R:", matrixR);
    }

    printf("\n\nC is Strassen result and R is ijk_matmul result. They should be equal. ");

    if (!matrix_equal(matrixC, matrixR))
        printf("Unfortunately they are not.\n\n");
    else
        printf("And indeed they are.\n\n");

    return 0;
}
