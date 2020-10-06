#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#define TIME(t_i,t_f) ((double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0) - \
                      ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);
#define RUNS 100
#define MATRIX_SIZE 8096

void random_matrix(float **M) {
    for (unsigned int i = 0; i < MATRIX_SIZE; ++i)
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) { 
            M[i][j] = (float)rand() / (float)RAND_MAX;
        }
}

int row_sum(float **m) {
    float result = 0;
    for (unsigned int i = 0; i < MATRIX_SIZE; ++i) 
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) {
            result += m[i][j];
        }
    
    return result;
} 

int column_sum(float **m) {
    float result = 0;
    for (unsigned int i = 0; i < MATRIX_SIZE; ++i) 
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) {
            result += m[j][i];
        }
    
    return result;
}

void benchmark(double *row_time, double *col_time) {
    // Generate random matrix
    float ** matrix; 
    matrix = (float**) malloc(MATRIX_SIZE*sizeof(float*));
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        matrix[i] = (float*) malloc(MATRIX_SIZE*sizeof(float)); 
    }
    random_matrix(matrix);

    struct timeval t_i, t_f;

    // Evaluate row_sum timing
    gettimeofday(&t_i, NULL);
    int res_row = row_sum(matrix);
    gettimeofday(&t_f, NULL);
    double t_sgetrf_row = TIME(t_i,t_f);

    // Evaluate column_sum timing
    gettimeofday(&t_i, NULL);
    int res_col = column_sum(matrix); 
    gettimeofday(&t_f, NULL);
    double t_sgetrf_col = TIME(t_i,t_f);

    assert(res_col == res_row);

    *row_time += t_sgetrf_row;
    *col_time += t_sgetrf_col;

    // Liberar memoria
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}


int main(int argc, char *argv[]) {
    srand(42); // Set random seed
    double row_time = 0;
    double col_time = 0;

    for (unsigned int j = 0; j < RUNS; j++) {
        benchmark(&row_time, &col_time);
    }
    printf("Results for iterating on a %ix%i matrix:\n\n", MATRIX_SIZE, MATRIX_SIZE);
    printf("Iterating through rows: %f ms\n", row_time/RUNS);
    printf("Iterating through cols: %f ms\n", col_time/RUNS);

    return 0;
}
