#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#define TIME(t_i,t_f) ((double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0) - \
                      ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);
#define RUNS 2
#define MATRIX_SIZE 4096

void random_matrix(float **M) {
    for (unsigned int i = 0; i < MATRIX_SIZE; ++i)
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) { 
            M[i][j] = (float)rand() / (float)RAND_MAX;
        }
}

void zero_matrix(float **M) {
    for (unsigned int i = 0; i < MATRIX_SIZE; ++i)
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) { 
            M[i][j] = (float)rand() / (float)RAND_MAX;
        }
}

void normal_mult(float **A, float **B, float **C, int m, int p, int n) {
    // Multiplying A and B and storing in C.
    for (unsigned int i = 0; i < m; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
            for (unsigned int k = 0; k < p; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void row_mult(float **A, float **B, float **C, int m, int p, int n) {
    // Multiplying A and B and storing in C.
    for (unsigned int i = 0; i < m; ++i) {
        for (unsigned int k = 0; k < p; ++k) {
            for (unsigned int j = 0; j < n; ++j) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

/* Ej C)
*   Construir la multiplicacion de matrices "por bloques" de tamanio variable (nb = 64, 128, 256, 512)
*   Las matrices no tienen porque ser cuadradas.
*   C = A x B. C_{ij} = Sum{k=1}{p}(A_{ik} x B_{kj})
*/
void mult_por_bloques(float **A, float **B, float **C, int m, int p, int n, int nb) {
    int blockSize = nb;

    // Initializing elements of matrix mult to 0.
    for (unsigned int i = 0; i < m; ++i)
        for (unsigned int j = 0; j < n; ++j)
            C[i][j] = 0;

    for (unsigned int bi=0; bi<m; bi+=blockSize)
        for (unsigned int bj=0; bj<n; bj+=blockSize)
            for (unsigned int bk=0; bk<p; bk+=blockSize)
                for (unsigned int i=0; i<blockSize; i++)
                    for (unsigned int k=0; k<blockSize; k++)
                        for (unsigned int j=0; j<blockSize; j++)
                            C[bi+i][bj+j] += A[bi+i][bk+k]*B[bk+k][bj+j];
}

void initialize_data(float ***A, float ***B, float ***C) {
    *A = (float**) malloc(MATRIX_SIZE*sizeof(float*));
    for (int i = 0; i < MATRIX_SIZE; ++i) (*A)[i] = (float*) malloc(MATRIX_SIZE*sizeof(float)); 

    *B = (float**) malloc(MATRIX_SIZE*sizeof(float*));
    for (int i = 0; i < MATRIX_SIZE; ++i) (*B)[i] = (float*) malloc(MATRIX_SIZE*sizeof(float)); 

    *C = (float**) malloc(MATRIX_SIZE*sizeof(float*));
    for (int i = 0; i < MATRIX_SIZE; ++i) (*C)[i] = (float*) malloc(MATRIX_SIZE*sizeof(float)); 

    random_matrix(*A);
    random_matrix(*B);
    zero_matrix(*C);
}


void free_data(float ***A, float ***B, float ***C) {
    for (int i = 0; i < MATRIX_SIZE; ++i) free((*A)[i]);
    for (int i = 0; i < MATRIX_SIZE; ++i) free((*B)[i]);
    for (int i = 0; i < MATRIX_SIZE; ++i) free((*C)[i]);
    free(*A);
    free(*B);
    free(*C);
}


void benchmark(double *row_time, double *col_time) {
    // Generate random matrix
    float **A, **B, **C;
    initialize_data(&A, &B, &C);

    struct timeval t_i, t_f;

    // Evaluate row_sum timing
    gettimeofday(&t_i, NULL);
    normal_mult(A, B, C, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE);
    gettimeofday(&t_f, NULL);
    double t_sgetrf_row = TIME(t_i,t_f);

    // Evaluate column_sum timing
    gettimeofday(&t_i, NULL);
    row_mult(A, B, C, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE);
    gettimeofday(&t_f, NULL);
    double t_sgetrf_col = TIME(t_i,t_f);

    *row_time += t_sgetrf_row;
    *col_time += t_sgetrf_col;

    free_data(&A, &B, &C);
}


int main(int argc, char *argv[]) {
    srand(42); // Set random seed
    double row_time = 0;
    double col_time = 0;

    for (unsigned int j = 0; j < RUNS; j++) {
        benchmark(&row_time, &col_time);
    }
    printf("Results for multiplicating two %ix%i matrixes:\n\n", MATRIX_SIZE, MATRIX_SIZE);
    printf("Multiplication through rows/columns: %f ms\n", row_time/RUNS);
    printf("Multiplication through rows/rows: %f ms\n", col_time/RUNS);

    return 0;
}
