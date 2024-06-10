#include "kmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define ITER 200
#define ERROR_K "Invalid number of clusters!"
#define ERROR_N "Invalid number of points!"
#define ERROR_d "Invalid dimension of point!"
#define ERROR_iter "Invalid maximum iteration!"

struct Cluster
{
    double *centroid;
    int size;
};
typedef struct Cluster Cluster;

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

/*used to free memory*/
void free_clusters(Cluster *clusters, int k)
{
    int i;
    for (i = 0; i < k; ++i)
    {
        free(clusters[i].centroid);
    }
    free(clusters);
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

/*finds the closest cluster to a given vector by calculating the euc dist from it*/
int find_closest_centroid_index(double **centroids, double *v, int k, int d)
{
    double min_dist = euc_l2(v, centroids[0], d);
    int index = 0;
    int i;
    for (i = 0; i < k; i++)
    {
        double dist = euc_l2(v, centroids[i], d);
        if (dist < min_dist)
        {
            index = i;
            min_dist = dist;
        }
    }
    return index;
}

/*calculates the average of all vectors in a cluster*/
double *calc_centroid_average(Cluster cluster, int d)
{
    int i, size;
    double *centroid, *averaged_vector;
    centroid = cluster.centroid;
    size = cluster.size;
    averaged_vector = calloc(d, sizeof(double));
    if (size == 0)
    {
        return averaged_vector;
    }
    for (i = 0; i < d; i++)
    {
        averaged_vector[i] = (double)centroid[i] / size;
    }
    return averaged_vector;
}

/*checks if the last change of centroids is less than epsilon for each*/
int check_centroid_convergence(double **centroids, double **new_centroids, int k, int d, double epsilon)
{
    int i, convergent_centroids;
    double *centroid_old, *centroid_new;
    convergent_centroids = 0;
    for (i = 0; i < k; i++)
    {
        centroid_old = centroids[i];
        centroid_new = new_centroids[i];
        if (euc_l2(centroid_old, centroid_new, d) <= epsilon)
            convergent_centroids += 1;
    }
    return (convergent_centroids == k);
}

/*adds vector to a cluster*/
void add_vector_to_centroid(Cluster *clus, double const vec[], int d)
{
    int i;
    double updated_entry_i;
    for (i = 0; i < d; i++)
    {
        updated_entry_i = clus->centroid[i] + vec[i];
        clus->centroid[i] = updated_entry_i;
    }
    clus->size = clus->size + 1;
}

double **k_means(int k, int n, int d, double epsilon, int iter, int *centroids_initial_indexes, double **data)
{
    /*creates the initial centroids according to the first k vectors given in the input*/
    int i, j, data_i, convergence, closest_centroid_index;
    double **centroids, **updated_centroids;
    double *x;
    centroids = sub_matrix_k(data, k, d, centroids_initial_indexes);
    for (i = 0; i < iter; i++)
    {
        struct Cluster *new_centroids = (struct Cluster *)createArray(k, sizeof(Cluster));
        for (j = 0; j < k; j++)
        {
            new_centroids[j].size = 0;
            new_centroids[j].centroid = calloc(d, sizeof(double));
        }
        for (data_i = 0; data_i < n; data_i++)
        {
            x = data[data_i];
            closest_centroid_index = find_closest_centroid_index(centroids, x, k, d);
            add_vector_to_centroid(&new_centroids[closest_centroid_index], x, d);
        }
        updated_centroids = (double **)createArray(k, sizeof(double *));
        for (j = 0; j < k; j++)
        {
            updated_centroids[j] = calc_centroid_average(new_centroids[j], d);
        }
        convergence = check_centroid_convergence(updated_centroids, centroids, k, d, epsilon);
        free_matrix(centroids, k);
        free_clusters(new_centroids, k);
        centroids = updated_centroids;
        if (convergence)
            break;
    }
    return centroids;
}

int myFunction()
{
    return 52;
}

// int main(int argc, char *argv[])
// {
//     int iter, i, j, K, n, d;
//     char *end1, *end2, *end3, *end4;
//     double **output, **data;
//     if (isStringDigit((char *)argv[3]) == 0)
//     {
//         printf(ERROR_d);
//         return 1;
//     }
//     if (isStringDigit((char *)argv[1]) == 0)
//     {
//         printf(ERROR_K);
//         return 1;
//     }
//     if (isStringDigit((char *)argv[2]) == 0)
//     {
//         printf(ERROR_N);
//         return 1;
//     }
//     K = strtol(argv[1], &end1, 10);
//     n = strtol(argv[2], &end2, 10);
//     d = strtol(argv[3], &end3, 10);
//     data = (double **)createArray(n, sizeof(double *));
//     for (i = 0; i < n; i++)
//     {
//         data[i] = (double *)createArray(d, sizeof(double));
//     }
//     for (i = 0; i < n; i++)
//     {
//         for (j = 0; j < d; j++)
//         {
//             char c;
//             double m;
//             scanf("%lf%c", &m, &c);
//             data[i][j] = m;
//         }
//     }
//     if (argc >= 5)
//     {
//         if (isStringDigit((char *)argv[4]) == 0)
//         {
//             printf(ERROR_iter);
//             free_matrix(data, n);
//             return 1;
//         }
//         else
//         {
//             iter = strtol(argv[4], &end4, 10);
//         }
//     }
//     else
//     {
//         iter = ITER;
//     }
//     if (n <= 1)
//     {
//         printf(ERROR_N);
//         free_matrix(data, n);
//         return 1;
//     }
//     if (K <= 1 || n <= K)
//     {
//         printf(ERROR_K);
//         free_matrix(data, n);
//         return 1;
//     }
//     if (d < 1)
//     {
//         printf(ERROR_d);
//         free_matrix(data, n);
//         return 1;
//     }
//     if (iter <= 1 || iter >= 1000)
//     {
//         printf(ERROR_iter);
//         free_matrix(data, n);
//         return 1;
//     }

//     output = k_means(K, n, d, iter, data);

//     for (i = 0; i < K; i++)
//     {
//         for (j = 0; j < d; j++)
//         {
//             printf("%.4f", output[i][j]);
//             if (j < d - 1)
//                 printf(",");
//         }
//         printf("\n");
//     }
//     free_matrix(data, n);
//     free_matrix(output, K);

//     return 0;
// }