#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define ITER 200
#define BETA 0.5
#define EPSILON 1e-4
/*
Calculates the similarity matrix of X. 
Input: a matrix X with N vectors of size d represented by a 1d array of size N*d
Output: A matrix of dimension N*N represented by a 1d array of size N*N
*/
double* sym(double* X, int N, int d){
    double* A = (double*) createArray(N*N, sizeof(double));
    int i, j = 0;
    for (i; i<N; i++){
        for(j; j<N; j++){
            double a_i_j;
            if(i !=j){
                a_i_j = exp(-pow(euc_l2(X + d*i, X + d*j, d), 2)/2 );
            } else{
                a_i_j = 0;
            }
            X[get_1d_index(i,j,d)] = a_i_j;
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
    int i = 0;
    for(i; i<N; i++){
        double d_i = 0;
        int j = 0;
        for (j; j<N; j++)
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
    double* D_invsqrt = calculate_inverse_square_root(D, N);
    double* D_invsqrtA = createArray(N*N, sizeof(double));
    matrix_multiply(D_invsqrt, A, D_invsqrt, N, N, N);
    double* normalized = createArray(N*N, sizeof(double));
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

    for(iter; iter < ITER; iter ++){
        matrix_multiply(W, H_curr, WH, N, N, k);
        transpose_matrix(H_curr, H_transposed, N, k);
        matrix_multiply(H_curr, H_transposed, HH_T, N, k, N);
        matrix_multiply(HH_T, H_curr, HH_TH, N, N, k);
        int i = 0;
        int j = 0;
        for(i; i<N; i++){
            for(j; j<k; j++){
                double alpha = (1 - BETA  + BETA*( WH[get_1d_index(i, j, k)] / HH_TH[get_1d_index(i, j, k)]));
                H_next[i*k + j] = H_curr[i*k+j] * alpha;
            }
        }
        double diff = frobenius_norm(H_next, H_curr, N, k);
        H_curr = H_next;
        if (diff < EPSILON)
            break;
    }

    free(WH);
    free(H_transposed);
    free(HH_T);
    free(HH_TH);
    free(H_next);
    return H_curr;
}

void matrix_multiply(double* A, double* B, double* C, int N, int d, int K) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < K; j++) 
            for (int k = 0; k < d; k++)
                C[get_1d_index(i,j, K)] += A[get_1d_index(i,k, d)] * B[get_1d_index(k, j, K)];
}

void transpose_matrix(int* A, double* transposed, int N, int d) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < d; j++) {
            transposed[j * N + i] = A[i * d + j];
        }
    }
}



double* calculate_inverse_square_root(double* D, int N) {
    double* D_inv_sqrt = (double*) createArray(N*N, sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
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


double frobenius_norm(double* A, double* B, int N, int K) {
    double sum = 0.0;

    for (int i = 0; i < N * K; i++) {
        double diff = A[i] - B[i];
        sum += diff * diff;
    }

    return sum;
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
