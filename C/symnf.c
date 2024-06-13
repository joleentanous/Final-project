#include "kmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define ITER 200

/*
Calculates the similarity matrix of X. 
Input: a matrix X with N vectors of size d represented by a 1d array of size N*d
Output: A matrix of dimension N*N represented by a 1d array of size N*N
*/
double* sym(double* X, int N, int d){
    double* A = (double*) createArray(N*N, sizeof(double))
    int i,j = 0
    for (i; i<N; i++){
        for(j; j<N; j++){
            double a_i_j;
            if(i !=j){
                a_i_j = exp(-pow(euc_l2(X + d*i, X + d*j, d), 2)/2 )
            } else{
                a_i_j = 0
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
    double* D_invsqrtA = matrix_multiply(D_invsqrt, A, N);
    double* normalized = matrix_multiply(D_invsqrtA, D_invsqrt, N);
    free(D_invsqrt);
    free(D_invsqrtA);
    return normalized;

}

void matrix_multiply(double* A, double* B, int N) {
    double* C = (double*) createArray(N*N, sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) 
            for (int k = 0; k < N; k++) 
                C[get_1d_index(i,j,N)] += A[get_1d_index(i,k, N)] * B[get_1d_index(k, j, N)];
    return C;

}

void calculate_inverse_square_root(double* D, int N) {
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


/*check if a string represents an integer*/
int isStringDigit(const char *str)
{
    /* Iterate through each character in the string*/
    int i;
    for (i = 0; str[i] != '\0'; i++)
    {
        /* Check if the character is not a digit*/
        if (!isdigit(str[i]))
        {
            /* If any character is not a digit, return false*/
            return 0;
        }
    }
    /* If all characters are digits, return true*/
    return 1;
}


/*creates an array*/
void *createArray(int n, int size)
{
    void *array = calloc(n * size);
    if (array == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    return array;
}

/*creates a submatrix*/
double **sub_matrix_k(double **matrix, int k, int d, int *indexes)
{
    int i;
    int j;
    double **sub_array;
    sub_array = malloc(k * sizeof(double *));
    for (i = 0; i < k; i++)
    {
        sub_array[i] = (double *)createArray(d, sizeof(double));
        for (j = 0; j < d; j++)
        {
            sub_array[i][j] = matrix[indexes[i]][j];
        }
    }

    return sub_array;
}

/*used to free memory*/
void free_matrix(double **matrix, int k)
{
    int i;
    for (i = 0; i < k; i++)
        free(matrix[i]);
    free(matrix);
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

double* extract_vector(double* X, int index,  int d){
    double* vector = (double*) createArray(d, sizeof(double));
    int i = 0
    int start = index*d;
    for(i; i < d; i++){
        vector[i] = X[start+i];
    }
    return vector;
}

int get_1d_index(int i, int j, int d){
    return i*d+j;
}