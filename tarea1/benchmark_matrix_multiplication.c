#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#define TIME(t_i,t_f) ((double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0) - \
                      ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);
#define RUNS 1
#define MATRIX_SIZE 1024

void random_matrix(float **M) {
    for (unsigned int i = 0; i < MATRIX_SIZE; ++i)
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) { 
            M[i][j] = (float)rand() / (float)RAND_MAX;
        }
}

void zero_matrix(float **M, int size) {
    for (unsigned int i = 0; i < size; ++i)
        for (unsigned int j = 0; j < size; ++j) { 
            M[i][j] = 0;
        }
}

void initialize_zero_matrix(float ***C, int size)  {
    *C = (float**) malloc(size*sizeof(float*));
    for (int i = 0; i < size; ++i) (*C)[i] = (float*) malloc(size*sizeof(float));
}

void initialize_random_data(float ***A, float ***B) {
    *A = (float**) malloc(MATRIX_SIZE*sizeof(float*));
    for (int i = 0; i < MATRIX_SIZE; ++i) (*A)[i] = (float*) malloc(MATRIX_SIZE*sizeof(float)); 

    *B = (float**) malloc(MATRIX_SIZE*sizeof(float*));
    for (int i = 0; i < MATRIX_SIZE; ++i) (*B)[i] = (float*) malloc(MATRIX_SIZE*sizeof(float));

    random_matrix(*A);
    random_matrix(*B);
}

void free_matrix(float ***A, int size) {
    for (int i = 0; i < size; ++i) free((*A)[i]);
    free(*A);
}

float** add(float** A, float** B, int n) {
    float** temp;
    initialize_zero_matrix(&temp, n);
    for(unsigned int i=0; i<n; ++i)
        for(unsigned int j=0; j<n; ++j)
            temp[i][j] = A[i][j] + B[i][j];
    return temp;
}

float** subtract(float** A, float** B, int n) {
    float** temp;
    initialize_zero_matrix(&temp, n);
    for(unsigned int i=0; i<n; i++)
        for(unsigned int j=0; j<n; j++)
            temp[i][j] = A[i][j] - B[i][j];
    return temp;
}

// Complexity O^3
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

// Complexity O^3
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

// Complexity O^3
void col_mult(float **A, float **B, float **C, int m, int p, int n) {
    // Multiplying A and B and storing in C.
    for (unsigned int j = 0; j < n; ++j) {
        for (unsigned int k = 0; k < p; ++k) {
            for (unsigned int i = 0; i < m; ++i) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


// Complexity O^2.808
float** strassen_mult(float **A, float **B, int n) {
    float** C;
    initialize_zero_matrix(&C, n);
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return C;
    }

    int k = n/2;
    float **A11, **A12, **A21, **A22, **B11, **B12, **B21, **B22;
    initialize_zero_matrix(&A11, k);
    initialize_zero_matrix(&A12, k);
    initialize_zero_matrix(&A21, k);
    initialize_zero_matrix(&A22, k);
    initialize_zero_matrix(&B11, k);
    initialize_zero_matrix(&B12, k);
    initialize_zero_matrix(&B21, k);
    initialize_zero_matrix(&B22, k);

    for(unsigned int i=0; i<k; ++i)
        for(unsigned int j=0; j<k; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][k+j];
            A21[i][j] = A[k+i][j];
            A22[i][j] = A[k+i][k+j];
            B11[i][j] = B[i][j];
            B12[i][j] = B[i][k+j];
            B21[i][j] = B[k+i][j];
            B22[i][j] = B[k+i][k+j];
        }

    float** P1 = strassen_mult(A11, subtract(B12, B22, k), k);
    float** P2 = strassen_mult(add(A11, A12, k), B22, k);
    float** P3 = strassen_mult(add(A21, A22, k), B11, k);
    float** P4 = strassen_mult(A22, subtract(B21, B11, k), k);
    float** P5 = strassen_mult(add(A11, A22, k), add(B11, B22, k), k);
    float** P6 = strassen_mult(subtract(A12, A22, k), add(B21, B22, k), k);
    free_matrix(&A12, k);
    free_matrix(&A22, k);
    free_matrix(&B21, k);
    free_matrix(&B22, k);

    float** P7 = strassen_mult(subtract(A11, A21, k), add(B11, B12, k), k);
    free_matrix(&A11, k);
    free_matrix(&A21, k);
    free_matrix(&B11, k);
    free_matrix(&B12, k);

    float** C11 = subtract(add(add(P5, P4, k), P6, k), P2, k);
    free_matrix(&P6, k);

    float** C12 = add(P1, P2, k);
    free_matrix(&P2, k);

    float** C21 = add(P3, P4, k);
    free_matrix(&P4, k);

    float** C22 = subtract(subtract(add(P5, P1, k), P3, k), P7, k);
    free_matrix(&P1, k);
    free_matrix(&P3, k);
    free_matrix(&P5, k);
    free_matrix(&P7, k);

    for(unsigned int i=0; i<k; ++i)
        for(unsigned int j=0; j<k; ++j) {
            C[i][j] = C11[i][j];
            C[i][j+k] = C12[i][j];
            C[k+i][j] = C21[i][j];
            C[k+i][k+j] = C22[i][j];
        }
    
    free_matrix(&C11, k);
    free_matrix(&C12, k);
    free_matrix(&C21, k);
    free_matrix(&C22, k);

    return C;
}

void assert_matrix_equality(float **A, float **B) {
    for (unsigned int i = 0; i < MATRIX_SIZE; ++i)
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j)
            if (A[i][j] != B[i][j])
                printf("%f, %f\n", A[i][j], B[i][j]); //assert(A[i][j] == B[i][j]);
}

void benchmark(double *normal_time, double *row_time, double *col_time, double *strassen_time) {
    // Generate random matrix
    float **A, **B, **C;
    initialize_random_data(&A, &B);
    initialize_zero_matrix(&C, MATRIX_SIZE);

    struct timeval t_i, t_f;

    // Evaluate normal_mult timing
    gettimeofday(&t_i, NULL);
    normal_mult(A, B, C, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE);
    gettimeofday(&t_f, NULL);
    *normal_time += TIME(t_i,t_f);
    zero_matrix(C, MATRIX_SIZE);

    // Evaluate column_mult timing
    gettimeofday(&t_i, NULL);
    row_mult(A, B, C, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE);
    gettimeofday(&t_f, NULL);
    *row_time += TIME(t_i,t_f);
    zero_matrix(C, MATRIX_SIZE);

    // Evaluate column_mult timing
    gettimeofday(&t_i, NULL);
    col_mult(A, B, C, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE);
    gettimeofday(&t_f, NULL);
    *col_time += TIME(t_i,t_f);
    // zero_matrix(C, MATRIX_SIZE);

    // Evaluate strassen_mult timing
    gettimeofday(&t_i, NULL);
    float ** D = strassen_mult(A, B, MATRIX_SIZE);
    gettimeofday(&t_f, NULL);
    *strassen_time += TIME(t_i,t_f);
    // assert_matrix_equality(C, D);

    free_matrix(&A, MATRIX_SIZE);
    free_matrix(&B, MATRIX_SIZE);
    free_matrix(&C, MATRIX_SIZE);
}


int main(int argc, char *argv[]) {
    srand(42); // Set random seed
    double normal_time = 0;
    double row_time = 0;
    double col_time = 0;
    double strassen_time = 0;
    
    for (unsigned int j = 0; j < RUNS; j++) {
        benchmark(&normal_time, &row_time, &col_time, &strassen_time);
    }
    printf("Results for multiplicating two %ix%i matrixes:\n\n", MATRIX_SIZE, MATRIX_SIZE);
    printf("Multiplication through cols/cols: %f ms\n", col_time/RUNS);
    printf("Multiplication through rows/columns: %f ms\n", normal_time/RUNS);
    printf("Multiplication through rows/rows: %f ms\n", row_time/RUNS);
    printf("Multiplication using strassen algorithm: %f ms\n", strassen_time/RUNS);

    return 0;
}
