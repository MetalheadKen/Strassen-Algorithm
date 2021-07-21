
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "Matrix.c"

#define N 64

static inline Matrix * strassen(Matrix * dest, const Matrix * srcA, const Matrix * srcB, int length)
{
 
    if (length == 2) return Matrix_Multiply(dest, srcA, srcB);

    int len = length / 2;

     matrix_autofree Matrix  *  a11= Matrix_Initializer( len, len),  * a12= Matrix_Initializer( len, len),  * a21= Matrix_Initializer( len, len),  * a22= Matrix_Initializer( len, len),  
     * b11= Matrix_Initializer( len, len),  * b12= Matrix_Initializer( len, len),  * b21= Matrix_Initializer( len, len),  * b22= Matrix_Initializer( len, len),  * c11= Matrix_Initializer( len, len),  * c12= Matrix_Initializer( len, len),  * c21= Matrix_Initializer( len, len),  * c22= Matrix_Initializer( len, len),  *              
                    m1= Matrix_Initializer( len, len),  * m2= Matrix_Initializer( len, len),  * m3= Matrix_Initializer( len, len),  * m4= Matrix_Initializer( len, len),  * m5= Matrix_Initializer( len, len),  * m6= Matrix_Initializer( len, len),  * m7= Matrix_Initializer( len, len),  * temp1= Matrix_Initializer( len, len),  * temp2 = Matrix_Initializer( len, len) ;

    /* Divide matrix into four parts */
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            a11->values[i][j] = srcA->values[i][j];
            a12->values[i][j] = srcA->values[i][j + len];
            a21->values[i][j] = srcA->values[i + len][j];
            a22->values[i][j] = srcA->values[i + len][j + len];
            
            b11->values[i][j] = srcB->values[i][j];
            b12->values[i][j] = srcB->values[i][j + len];
            b21->values[i][j] = srcB->values[i + len][j];
            b22->values[i][j] = srcB->values[i + len][j + len];
        }
    }

    /* Calculate seven formulas of strassen Algorithm */
    strassen(m1, Matrix_Addition(temp1, a11, a22), Matrix_Addition(temp2, b11, b22), len);
    strassen(m2, Matrix_Addition(temp1, a21, a22), b11, len);
    strassen(m3, a11, Matrix_Subtract(temp1, b12, b22), len);
    strassen(m4, a22, Matrix_Subtract(temp1, b21, b11), len);
    strassen(m5, Matrix_Addition(temp1, a11, a12), b22, len);
    strassen(m6, Matrix_Subtract(temp1, a21, a11), Matrix_Addition(temp2, b11, b12), len);
    strassen(m7, Matrix_Subtract(temp1, a12, a22), Matrix_Addition(temp2, b21, b22), len);

    /* Merge the answer of matrix dest */
    /* c11 = m1 + m4 - m5 + m7 = m1 + m4 - (m5 - m7) */
    Matrix_Subtract(c11, Matrix_Addition(temp1, m1, m4), Matrix_Subtract(temp2, m5, m7));
    Matrix_Addition(c12, m3, m5);
    Matrix_Addition(c21, m2, m4);
    Matrix_Addition(c22, Matrix_Subtract(temp1, m1, m2), Matrix_Addition(temp2, m3, m6));
    
    /* Store the answer of matrix multiplication */
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            dest->values[i][j]              = c11->values[i][j];
            dest->values[i][j + len]        = c12->values[i][j];
            dest->values[i + len][j]        = c21->values[i][j];
            dest->values[i + len][j + len]  = c22->values[i][j];
        }
    }

    return dest;
}

static void matrix_print ( const char * prompt, Matrix * mx ) 
{
    assert( prompt && mx );
    printf("%s\n", prompt );
    for(int i = 0; i < mx->row  ; i++) {
        for (int j = 0; j < mx->column; j++) {
            printf(" %5d ", mx->values[i][j]);
        }
        printf("\n");
    }
}

static void matrix_foreach ( Matrix * mx, Matrix * (*callback)(Matrix *, unsigned, unsigned) ) 
{
    for(int i = 0; i < mx->row  ; i++) {
        for (int j = 0; j < mx->column; j++) {
            callback(mx,i,j);
        }
    }
}

Matrix * ordinal_as_val (Matrix * mx, unsigned i, unsigned j) 
{
    mx->values[i][j] = (i * mx->column + j) ;
    return mx ;
}

/*
--------------------------------------------------------------------------------------------------------------------
*/
int main(int argc, char *argv[])
{
    const unsigned matrix_length = N ;
  
    /* Check if dimensions of matrix is the power of two or not */
    if (ceil(log2(matrix_length)) != floor(log2(matrix_length)))
    {
        printf("\nERROR: square matrix side must be a power of 2. And %d is not.",matrix_length );
        return 0;
    }

    matrix_autofree Matrix * matrixA = Matrix_Initializer( matrix_length, matrix_length);
    matrix_autofree Matrix * matrixB = Matrix_Initializer( matrix_length, matrix_length);
    matrix_autofree Matrix * matrixC = Matrix_Initializer( matrix_length, matrix_length);
    
    matrix_foreach(matrixA, ordinal_as_val );
    matrix_foreach(matrixB, ordinal_as_val );

    /* Matrix multiplication */
    matrixC = strassen(matrixC, matrixA, matrixB, matrix_length);

matrix_print("Matrix A:", matrixA );
matrix_print("Matrix B:", matrixB );
matrix_print("Matrix C:", matrixC );

    return 0;
}
