#include <Python.h>
#define PY_SSIZE_T_CLEAN
#include "kmeans.h"
#include <stdio.h>
#include <math.h>

double *spread_matrix(double **matrix, int N, int d)
{
    double *arr = (double *)malloc(N * d * sizeof(double));
    if (arr == NULL)
    {
        printf("Memory allocation failed. Exiting.\n");
        return NULL;
    }

    int index = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < d; j++)
        {
            arr[index++] = matrix[i][j];
        }
    }

    return arr;
}

void printArray(int *arr, int length)
{
    printf("[");
    for (int i = 0; i < length; ++i)
    {
        printf("%d", arr[i]);
        if (i < length - 1)
        {
            printf(", ");
        }
    }
    printf("]\n");
}

void print2DArray(double **arr, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%f\t", arr[i][j]);
        }
        printf("\n");
    }
}

double *k_means_wrapper(int K, int N, int d, int iter, double epsilon, int *indices, double **data)
{
    double **centroids = k_means(K, N, d, epsilon, iter, indices, data);
    double *spread_centroids = spread_matrix(centroids, K, d);
    free_matrix(centroids, K);
    return spread_centroids;
}

void print_matrix(int **matrix, int rows, int cols)
{
    printf("Received matrix:\n");
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

PyObject *convert_array_to_python_list(double *arr, int N)
{
    PyObject *py_list = PyList_New(N);
    if (py_list == NULL)
    {
        printf("Error creating Python list. Exiting.\n");
        return NULL;
    }

    for (int i = 0; i < N; i++)
    {
        PyObject *py_double = Py_BuildValue("d", arr[i]);
        if (py_double == NULL)
        {
            printf("Error creating Python double object. Exiting.\n");
            Py_DECREF(py_list);
            return NULL;
        }
        PyList_SetItem(py_list, i, py_double);
    }

    return py_list;
}

double **convert_array_to_matrix(double *arr, int N, int d)
{
    int i;
    double **matrix = (double **)malloc(N * sizeof(double *));
    if (matrix == NULL)
    {
        printf("Memory allocation failed. Exiting.\n");
        return NULL;
    }
    for (i = 0; i < N; i++)
    {
        matrix[i] = (double *)malloc(d * sizeof(double));
    }
    for (i = 0; i < N * d; i++)
    {
        double item = arr[i];
        matrix[i / d][i % d] = item;
    }
    free(arr);
    return matrix;
}

void *convert_list_to_int_array(PyObject *list)
{
    PyObject *item;
    int n = PyObject_Length(list);
    if (n < 0)
    {
        return NULL;
    }
    // printf("%d\n", PyObject_Length(list));
    int *arr = (int *)malloc(n * sizeof(int));
    if (arr == NULL)
    {
        printf("Memory allocation failed. Exiting.\n");
        return NULL;
    }
    int i;
    for (i = 0; i < n; i++)
    {
        item = PyList_GetItem(list, i);
        arr[i] = PyLong_AsLong(item);
    }
    return arr;
}

void *convert_list_to_double_array(PyObject *list)
{
    PyObject *item;
    int n = PyObject_Length(list);
    if (n < 0)
    {
        return NULL;
    }
    // printf("%d\n", PyObject_Length(list));
    double *arr = (double *)malloc(n * sizeof(double));
    if (arr == NULL)
    {
        printf("Memory allocation failed. Exiting.\n");
        return NULL;
    }
    int i;
    for (i = 0; i < n; i++)
    {
        item = PyList_GetItem(list, i);
        arr[i] = PyFloat_AsDouble(item);
    }
    return arr;
}

void parse_input(PyObject *args, int *K, int *N, int *d, int *iter, double *epsilon, int **indexes, double **data_array)
{
    // PyObject *lst;
    PyObject *indexes_list_ptr;
    PyObject *data_list_ptr;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if (!PyArg_ParseTuple(args, "iiiidOO", K, N, d, iter, epsilon, &indexes_list_ptr, &data_list_ptr))
    {
        return NULL;
    }
    // *data = process_matrix_interface(data_list_ptr);
    *indexes = convert_list_to_int_array(indexes_list_ptr);
    *data_array = convert_list_to_double_array(data_list_ptr);
}

static PyObject *fit(PyObject *self, PyObject *args)
{
    int K, N, d, iter;
    double epsilon;
    int *initial_centroid_indexes;
    double *data_array;
    double **data_matrix;

    parse_input(args, &K, &N, &d, &iter, &epsilon, &initial_centroid_indexes, &data_array);
    data_matrix = convert_array_to_matrix(data_array, N, d);
    double *spread_centroids = k_means_wrapper(K, N, d, iter, epsilon, initial_centroid_indexes, data_matrix);
    return convert_array_to_python_list(spread_centroids, K * d);
}

static PyMethodDef geoMethods[] = {
    {"fit",            /* the Python method name that will be used */
     (PyCFunction)fit, /* the C-function that implements the Python function and returns static PyObject*  */
     METH_VARARGS,     /* flags indicating parameters
accepted for this function */
     PyDoc_STR("Calculates k_means")},
    /*  The docstring for the function */
    {NULL, NULL, 0, NULL}, /* The last entry must be all NULL as shown to act as a
                              sentinel. Python looks for this entry to know that all
                                of the functions for the module have been defined. */
};

static struct PyModuleDef geomodule = {
    PyModuleDef_HEAD_INIT,
    "kmeansmodule", /* name of module */
    NULL,           /* module documentation, may be NULL */
    -1,             /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    geoMethods      /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_kmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&geomodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}