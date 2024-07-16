#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define ITER 200
#define BETA 0.5
#define EPSILON 1e-4

double* sym(double* X, int N, int d);
void printMatrix(double* matrix, int rows, int cols);
void *createArray(int n, int size);
int get_1d_index(int i, int j, int d);
double euc_l2(double *v1, double *v2, int d);
double* calculate_inverse_square_root(double* D, int N);
void matrix_multiply(double* A, double* B, double* C, int N, int d, int K);
double* ddg(double* A, int N);
double* norm(double* A, double* D, int N);
void transpose_matrix(double* A, double* transposed, int N, int d);
double frobenius_norm(double* A, double* B, int N, int K);
double* symnmf(double* W, double* H, int N, int k);
void copyArray(double* A, double* B, int size);

int main() {
    /* 
    int N = 3;
    int d = 2;

    double X[] = {
        1.0, 2.0,
        3.0, 4.0,
        5.0, 6.0
    };

    double* similarityMatrix;
    double* ddg_matrix;
    double* normalized_matrix;

    similarityMatrix = sym(X, N, d);

    printf("Similarity Matrix:\n");
    printMatrix(similarityMatrix, N, N);

    ddg_matrix = ddg(similarityMatrix, N);

    printf("DDG Matrix:\n");
    printMatrix(ddg_matrix, N, N);

    printf("Normalized Matrix:\n");
    normalized_matrix = norm(similarityMatrix, ddg_matrix, N);
    printMatrix(normalized_matrix, N, N);
    free(similarityMatrix);
    free(ddg_matrix);
    */

    int N = 3; 
    int k = 2;

    double W[] = {
        0.000000, 0.707105, 0.000006,
        0.707105, 0.000000, 0.707105,
        0.000006, 0.707105, 0.000000
    };

    double H[] = {
        0.5, 0.6,
        0.3, 0.8,
        0.9, 0.1
    };

    double* optimizedH = symnmf(W, H, N, k);

    printf("Optimized H Matrix:\n");
    printMatrix(optimizedH, N, k);

    /*free(optimizedH);*/

    return 0;
}

/*
Calculates the similarity matrix of X. 
Input: a matrix X with N vectors of size d represented by a 1d array of size N*d
Output: A matrix of dimension N*N represented by a 1d array of size N*N
*/
double* sym(double* X, int N, int d){
    double* A = (double*) createArray(N*N, sizeof(double));
    int i,j;
    for (i = 0; i<N; i++){
        double a_i_j;
        for(j=0; j<N; j++){
            if(i !=j){
                a_i_j = exp(-pow(euc_l2(X + d*i, X + d*j, d), 2)/2 );
            } else{
                a_i_j = 0;
            }
            A[get_1d_index(i,j,N)] = a_i_j;
        }
    }
    return A;
}


/*
Calculates the diagonal degree Matrix
Input: 1d array representing the similarity matrix A, the dimension N
Output: 1d array representing the diagonal degree Matrix D based on A
*/
double* ddg(double* A, int N){
    double* D = (double*) createArray(N*N, sizeof(double));
    int i;
    for(i = 0; i<N; i++){
        double d_i = 0;
        int j;
        for (j=0; j<N; j++)
            d_i += A[get_1d_index(i, j, N)];
        D[get_1d_index(i, i, N)] = d_i;
    }
    return D;
}



/*
Calculates the  normalized similarity matrix
Input: two 1d arrays representing the similarity matrix A and the ddg matrix D , the dimension N
Output: 1d array representing the  normalized similarity matrix D^-1/2AD^-1/2
*/
double* norm(double* A, double* D, int N){
    double* D_invsqrt;
    double* D_invsqrtA;
    double* normalized;
    D_invsqrt = calculate_inverse_square_root(D, N);
    D_invsqrtA = createArray(N*N, sizeof(double));
    matrix_multiply(D_invsqrt, A, D_invsqrtA, N, N, N);
    normalized = createArray(N*N, sizeof(double));
    matrix_multiply(D_invsqrtA, D_invsqrt, normalized, N, N, N);
    free(D_invsqrt);
    free(D_invsqrtA);
    return normalized;
}

/*
Optimizing the randomly generated H
Input: W - 1d array representation of the NxN normalized similarity matrix 
H - 1d array representation of a Nxk randomly initialized matrix
Output: Optimized H
*/
double* symnmf(double* W, double* H, int N, int k ){
    int iter = 0;
    double* H_curr = H;
    double* H_next = createArray(N*k, sizeof(double));
    double* WH = createArray(N*k, sizeof(double));
    double* H_transposed = createArray(N*k, sizeof(double));
    double* HH_T = createArray(N*N, sizeof(double));
    double* HH_TH = createArray(N*k, sizeof(double));

    for(iter=0; iter < ITER; iter ++){
        int i;
        int j;
        double diff;
        matrix_multiply(W, H_curr, WH, N, N, k);
        transpose_matrix(H_curr, H_transposed, N, k);
        matrix_multiply(H_curr, H_transposed, HH_T, N, k, N);
        matrix_multiply(HH_T, H_curr, HH_TH, N, N, k);
        for(i=0; i<N; i++){
            for(j=0; j<k; j++){
                double alpha = (1 - BETA  + BETA*( WH[get_1d_index(i, j, k)] / HH_TH[get_1d_index(i, j, k)]));
                H_next[i*k + j] = H_curr[i*k+j] * alpha;
            }
        }
        diff = frobenius_norm(H_next, H_curr, N, k);
        copyArray(H_next, H_curr, N*k);
        if (diff < EPSILON)
            break;
    }

    free(WH);
    free(H_transposed);
    free(HH_T);
    free(HH_TH);
    /*free(H_next);*/
    return H_curr;
}

void matrix_multiply(double* A, double* B, double* C, int N, int d, int K) {
    int i,j,k;
    for (i = 0; i < N; i++)
        for (j = 0; j < K; j++) 
            for (k = 0; k < d; k++)
                C[get_1d_index(i,j, K)] += A[get_1d_index(i,k, d)] * B[get_1d_index(k, j, K)];
}

double* calculate_inverse_square_root(double* D, int N) {
    double* D_inv_sqrt = (double*) createArray(N*N, sizeof(double));
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            int i_j_index = get_1d_index(i, j, N);
            if (i == j) {
                D_inv_sqrt[i_j_index] = 1.0 / sqrt(D[i_j_index]);
            } else {
                D_inv_sqrt[i_j_index] = 0.0;
            }
        }
    }
    return D_inv_sqrt;
}



void transpose_matrix(double* A, double* transposed, int N, int d) {
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++) {
            transposed[j * N + i] = A[i * d + j];
        }
    }
}

double frobenius_norm(double* A, double* B, int N, int K) {
    double sum = 0.0;
    int i;
    for (i = 0; i < N * K; i++) {
        double diff = A[i] - B[i];
        sum += diff * diff;
    }

    return sum;
}

void copyArray(double* A, double* B, int size) {
    int i;
    for (i = 0; i < size; i++) {
        B[i] = A[i];
    }
}


/*creates an array*/
void *createArray(int n, int size)
{
    void *array = calloc(n, size);
    if (array == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    return array;
}


/*calculates the euc distance between two given vectors*/
double euc_l2(double *v1, double *v2, int d)
{
    int i;
    double dist;
    dist = 0.0;
    for (i = 0; i < d; i++)
    {

        dist += pow(v1[i] - v2[i], 2.0);
    }
    return sqrt(dist);
}

int get_1d_index(int i, int j, int d){
    return i*d+j;
}


void printMatrix(double* matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            printf("%f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
}
