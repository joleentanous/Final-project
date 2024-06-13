#include "kmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define ITER 200

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
    void *array = malloc(n * size);
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