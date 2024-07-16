#include <Python.h>
#define PY_SSIZE_T_CLEAN
#include "symnmf.h"
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

void parse_input_symnmf(PyObject *args, double **W_array, double **H_array, int *N, int *k)
{
    
    PyObject *W_list_ptr;
    PyObject *H_list_ptr;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if (!PyArg_ParseTuple(args, "OOii", &W_list_ptr, &H_list_ptr, N, k))
    {
        return NULL;
    }
    
    *W_array = convert_list_to_double_array(W_list_ptr);
    *H_array = convert_list_to_double_array(H_list_ptr);
}


double *symnmf_wrapper(double **w, double **h, int N, int k)
{
    /*(W, H, N, k)*/
    double **symnmf = symnmf(W, H, N, k);
    double *symnmf_vectors = spread_matrix(symnmf, N, k);
    free_matrix(symnmf, N);
    return symnmf_vectors;
    
}

static PyObject *symnmf(PyObject *self, PyObject *args)
{
    int N, k;
    double *W_array, *H_array;
    double **W_matrix, **H_matrix;
    parse_input_norm(args, &W_array, &H_array, &N, &k);
    W_matrix = convert_array_to_matrix(W_array, N, N);
    H_matrix = convert_array_to_matrix(H_array, N, N);
    double *symnmf_vectors = symnmf_wrapper(W_matrix, H_matrix, N, k);
    return convert_array_to_python_list(symnmf_vectors, N * k);
}


void parse_input_sym(PyObject *args, double **data_array, int *N, int *d)
{
    PyObject *data_list_ptr;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if (!PyArg_ParseTuple(args, "Oii", &data_list_ptr, N, d))
    {
        return NULL;
    }
    *data_array = convert_list_to_double_array(data_list_ptr);
}

double *sym_wrapper(double **X, int N, int d)
{
    double **sym = sym(X, N, d);
    double *sym_vectors = spread_matrix(sym, N, N);
    free_matrix(sym, N);
    return sym_vectors;
    
}

static PyObject *sym(PyObject *self, PyObject *args)
{
    int N,d;
    double *X_array;
    double **X_matrix;
    parse_input_sym(args, &X_array, &N, &d);
    X_matrix = convert_array_to_matrix(X_array, N, d);
    double *sym_vectors = sym_wrapper(X_matrix, N, d);
    return convert_array_to_python_list(sym_vectors, N * N);
}


void parse_input_ddg(PyObject *args, double **data_array, int *N)
{
    PyObject *data_list_ptr;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if (!PyArg_ParseTuple(args, "Oi", &data_list_ptr, N))
    {
        return NULL;
    }
    *data_array = convert_list_to_double_array(data_list_ptr);
}


double *ddg_wrapper(double **A, int N)
{
    double **ddg = ddg(A, N);
    double *ddg_vectors = spread_matrix(ddg, N, N);
    free_matrix(ddg, N);
    return ddg_vectors;
    
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    int N;
    double *A_array;
    double **A_matrix;
    parse_input_ddg(args, &A_array, &N);
    A_matrix = convert_array_to_matrix(A_array, N, N);
    double *ddg_vectors = ddg_wrapper(A_matrix, N);
    return convert_array_to_python_list(ddg_vectors, N * N);
}



void parse_input_norm(PyObject *args, double **A_array, double **D_array, int *N)
{
    PyObject *A_list_ptr;
    PyObject *D_list_ptr;

    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if (!PyArg_ParseTuple(args, "OOi", &A_list_ptr, &D_list_ptr, N))
    {
        return NULL;
    }
    *A_array = convert_list_to_double_array(A_list_ptr);
    *D_array = convert_list_to_double_array(D_list_ptr);
}


double *norm_wrapper(double **A, double **A, int N)
{
    double **norm = norm(A, D, N);
    double *norm_vectors = spread_matrix(norm, N, N);
    free_matrix(norm, N);
    return norm_vectors;
    
}

static PyObject *norm(PyObject *self, PyObject *args)
{
    int N;
    double *A_array, *D_array;
    double **A_matrix, **D_matrix;
    parse_input_norm(args, &A_array, &D_array, &N);
    A_matrix = convert_array_to_matrix(A_array, N, N);
    D_matrix = convert_array_to_matrix(D_array, N, N);
    double *norm_vectors = norm_wrapper(A_matrix, D_matrix, N);
    return convert_array_to_python_list(norm_vectors, N * N);
    
}


static PyMethodDef geoMethods[] = {
    {"symnmf",            /* the Python method name that will be used */
     (PyCFunction)symnmf, /* the C-function that implements the Python function and returns static PyObject*  */
     METH_VARARGS,     /* flags indicating parameters accepted for this function */
     PyDoc_STR("Computes the whole process of symnmf")},
    /*  The docstring for the function */
    {"sym",
    (PyCFunction)sym,
    METH_VARARGS,
    PyDoc_STR("Computes sym part")},
    {"ddg",
    (PyCFunction)ddg,
    METH_VARARGS,
    PyDoc_STR("Computes ddg part")},
    {"norm",
    (PyCFunction)norm,
    METH_VARARGS,
    PyDoc_STR("Computes ddg part")},
    {NULL, NULL, 0, NULL}, /* The last entry must be all NULL as shown to act as a
                              sentinel. Python looks for this entry to know that all
                                of the functions for the module have been defined. */
};

static struct PyModuleDef geomodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule", /* name of module */
    NULL,           /* module documentation, may be NULL */
    -1,             /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    geoMethods      /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmfsmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}
