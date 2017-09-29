#ifdef _cplusplus
extern "C" {
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
}
#else
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
#endif

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

#define LOG2(LENGTH) \
    (8 * sizeof(unsigned int) - Count_Leading_Zero((LENGTH)) - 1)

typedef struct _Matrix {
    int row;
    int column;
    int **values;
    
    int **(*New)(int row, int column);
    void  (*Delete)(void *);
} Matrix;

int Count_Leading_Zero(unsigned int number);

int **Matrix_Allocate(int, int);
void  Matrix_Free(void *);

Matrix *Strassen(Matrix *, Matrix *, Matrix *, int);
Matrix *Matrix_Multiply(Matrix *, Matrix *, Matrix *);
Matrix *Matrix_Addition(Matrix *, Matrix *, Matrix *);
Matrix *Matrix_Subtract(Matrix *, Matrix *, Matrix *);

int main(int argc, char *argv[])
{
    int dimensions, matrix_length;

    printf("Please enter the dimensions of the matrix: ");
    scanf("%d", &dimensions);

    if (dimensions <= 0) { printf("The number you entered is invalid."); exit(1); }

    /* Check if dimensions of matrix is the power of two or not */
    if (dimensions & (dimensions - 1))
        matrix_length = 2 << LOG2((unsigned int) dimensions);
    else
        matrix_length = dimensions;

    Matrix matrixA = { .row = matrix_length, .column = matrix_length, .values = NULL, .New = Matrix_Allocate, .Delete = free };
    Matrix matrixB = { .row = matrix_length, .column = matrix_length, .values = NULL, .New = Matrix_Allocate, .Delete = free };
    Matrix matrixC = { .row = matrix_length, .column = matrix_length, .values = NULL, .New = Matrix_Allocate, .Delete = free };

    matrixA.values = matrixA.New(matrixA.row, matrixB.column);
    matrixB.values = matrixB.New(matrixB.row, matrixB.column);
    matrixC.values = matrixC.New(matrixC.row, matrixC.column);
    
    for (int i = 0; i < dimensions; i++) {
        for (int j = 0; j < dimensions; j++) {
            matrixA.values[i][j] = matrixB.values[i][j] = i * dimensions + j;
        }
    }

    /* Matrix multiplication */
    Strassen(&matrixC, &matrixA, &matrixB, matrixC.row);

    /* Print the answer of matrix multiplication */
    printf("Matrix A:\n");
    for(int i = 0; i < dimensions; i++) {
        for (int j = 0; j < dimensions; j++) {
            printf("%5d ", matrixA.values[i][j]);
        }
        printf("\n");
    }

    printf("\nMatrix B:\n");
    for(int i = 0; i < dimensions; i++) {
        for (int j = 0; j < dimensions; j++) {
            printf("%5d ", matrixB.values[i][j]);
        }
        printf("\n");
    }

    printf("\nMatrix C:\n");
    for(int i = 0; i < dimensions; i++) {
        for (int j = 0; j < dimensions; j++) {
            printf("%5d ", matrixC.values[i][j]);
        }
        printf("\n");
    }

    matrixA.Delete(matrixA.values);
    matrixB.Delete(matrixB.values);
    matrixC.Delete(matrixC.values);

    return 0;
}

int Count_Leading_Zero(unsigned int number)
{
    if (number == 0) return 32;
    
    int count = 0;
    if (number <= 0x0000FFFF) { count += 16; number <<= 16; }
    if (number <= 0x00FFFFFF) { count +=  8; number <<=  8; }
    if (number <= 0x0FFFFFFF) { count +=  4; number <<=  4; }
    if (number <= 0x3FFFFFFF) { count +=  2; number <<=  2; }
    if (number <= 0x7FFFFFFF) { count +=  1; number <<=  1; }

    return count;
}

Matrix *Strassen(Matrix *dest, Matrix *srcA, Matrix *srcB, int length)
{
    if (length == 2) return Matrix_Multiply(dest, srcA, srcB);

    int len = length / 2;

    Matrix a11, a12, a21, a22,              /* Matrix srcA divide four part */
           b11, b12, b21, b22,              /* Matrix srcB divide four part */
           c11, c12, c21, c22,              /* Matrix dest divide four part */
           m1, m2, m3, m4, m5, m6, m7,
           temp1, temp2;

    /* Initializer the matrix */
    MATRIX_INITIALIZER(a11, len, len);
    MATRIX_INITIALIZER(a12, len, len);
    MATRIX_INITIALIZER(a21, len, len);
    MATRIX_INITIALIZER(a22, len, len);
    MATRIX_INITIALIZER(b11, len, len);
    MATRIX_INITIALIZER(b12, len, len);
    MATRIX_INITIALIZER(b21, len, len);
    MATRIX_INITIALIZER(b22, len, len);
    MATRIX_INITIALIZER(c11, len, len);
    MATRIX_INITIALIZER(c12, len, len);
    MATRIX_INITIALIZER(c21, len, len);
    MATRIX_INITIALIZER(c22, len, len);
    MATRIX_INITIALIZER(m1, len, len);
    MATRIX_INITIALIZER(m2, len, len);
    MATRIX_INITIALIZER(m3, len, len);
    MATRIX_INITIALIZER(m4, len, len);
    MATRIX_INITIALIZER(m5, len, len);
    MATRIX_INITIALIZER(m6, len, len);
    MATRIX_INITIALIZER(m7, len, len);
    MATRIX_INITIALIZER(temp1, len, len);
    MATRIX_INITIALIZER(temp2, len, len);

    /* Divide matrix to four part */
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            a11.values[i][j] = srcA->values[i][j];
            a12.values[i][j] = srcA->values[i][j + len];
            a21.values[i][j] = srcA->values[i + len][j];
            a22.values[i][j] = srcA->values[i + len][j + len];
            
            b11.values[i][j] = srcB->values[i][j];
            b12.values[i][j] = srcB->values[i][j + len];
            b21.values[i][j] = srcB->values[i + len][j];
            b22.values[i][j] = srcB->values[i + len][j + len];
        }
    }

    /* Calculate seven formulas of Strassen Algorithm */
    Strassen(&m1, Matrix_Addition(&temp1, &a11, &a22), Matrix_Addition(&temp2, &b11, &b22), len);
    Strassen(&m2, Matrix_Addition(&temp1, &a21, &a22), &b11, len);
    Strassen(&m3, &a11, Matrix_Subtract(&temp1, &b12, &b22), len);
    Strassen(&m4, &a22, Matrix_Subtract(&temp1, &b21, &b11), len);
    Strassen(&m5, Matrix_Addition(&temp1, &a11, &a12), &b22, len);
    Strassen(&m6, Matrix_Subtract(&temp1, &a21, &a11), Matrix_Addition(&temp2, &b11, &b12), len);
    Strassen(&m7, Matrix_Subtract(&temp1, &a12, &a22), Matrix_Addition(&temp2, &b21, &b22), len);

    /* Merge the answer of matrix dest */
    /* c11 = m1 + m4 - m5 + m7 = m1 + m4 - (m5 - m7) */
    Matrix_Subtract(&c11, Matrix_Addition(&temp1, &m1, &m4), Matrix_Subtract(&temp2, &m5, &m7));
    Matrix_Addition(&c12, &m3, &m5);
    Matrix_Addition(&c21, &m2, &m4);
    Matrix_Addition(&c22, Matrix_Subtract(&temp1, &m1, &m2), Matrix_Addition(&temp2, &m3, &m6));
    
    /* Store the answer of matrix multiplication */
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            dest->values[i][j]              = c11.values[i][j];
            dest->values[i][j + len]        = c12.values[i][j];
            dest->values[i + len][j]        = c21.values[i][j];
            dest->values[i + len][j + len]  = c22.values[i][j];
        }
    }

    a11.Delete(a11.values); a12.Delete(a12.values); a21.Delete(a21.values); a22.Delete(a22.values);
    b11.Delete(b11.values); b12.Delete(b12.values); b21.Delete(b21.values); b22.Delete(b22.values);
    c11.Delete(c11.values); c12.Delete(c12.values); c21.Delete(c21.values); c22.Delete(c22.values);
    m1.Delete(m1.values);   m2.Delete(m2.values);   m3.Delete(m3.values);   m4.Delete(m4.values);
    m5.Delete(m5.values);   m6.Delete(m6.values);   m7.Delete(m7.values);
    temp1.Delete(temp1.values); temp2.Delete(temp2.values);

    return dest;
}

Matrix *Matrix_Multiply(Matrix *res, Matrix *a, Matrix *b)
{
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

Matrix *Matrix_Addition(Matrix *res, Matrix *a, Matrix *b)
{
    for (int i = 0; i < res->row; i++)
        for (int j = 0; j < res->column; j++)
            res->values[i][j] = a->values[i][j] + b->values[i][j];

    return res;
}

Matrix *Matrix_Subtract(Matrix *res, Matrix *a, Matrix * b)
{
    for (int i = 0; i < res->row; i++)
        for (int j = 0; j < res->column; j++)
            res->values[i][j] = a->values[i][j] - b->values[i][j];

    return res;
}

int **Matrix_Allocate(int row, int column)
{
    int **matrix = (int **) calloc(row + row * column, sizeof(**matrix));
    int *matrixTemp = (int *) (matrix + row);

    for (int i = 0; i < row; i++) {
        matrix[i] = matrixTemp;
        matrixTemp += column;
    }

    return matrix;
}
