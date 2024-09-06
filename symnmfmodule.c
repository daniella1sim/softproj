# define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "symnmf.h"


#define EPSILON 0.0001
#define MAX_ITERATIONS 300
#define BETA 0.5


void printCord(struct cord *c)
{
    while(c != NULL)
    {
        printf("%.4f, ", c->value);
        c = c->next;
    }
}


void printVector(struct vector *v)
{
    while(v != NULL)
    {
        printCord(v->cords);
        printf("\n");
        v = v->next;
    }
}


struct vector* PyObjectToLinkedList(PyObject *obj)
{
    PyObject *currObj;
    int rows, columns;

    struct vector *headVec, *currVec;
    struct cord *headCord, *currCord;

    headCord = (struct cord*)malloc(sizeof(struct cord));
    currCord = headCord;
    currCord->next = NULL;

    headVec = (struct vector*)malloc(sizeof(struct vector));
    currVec = headVec;
    currVec->next = NULL;

    rows = (int)PyList_Size(obj);
    if (rows == 0) columns = 0;
    else columns = (int) PyList_Size(PyList_GetItem(obj, 0));

    for (int i = 0; i < rows; i++)
    {
        currVec->cords = currCord;
        currObj = PyList_GetItem(obj, i);
        for(int j = 0; j < columns; j++)
        {
            currCord->value = PyFloat_AsDouble(PyList_GetItem(currObj, j));
            if (j == columns - 1) break;
            currCord = addNewCord(currCord);
        }
        if (i == rows - 1) break;
        currVec = addVector(currVec);
        headCord = createNewCord();
        currCord = headCord;
    }
    return headVec;
}


PyObject* matrixToPyObject(double** matrix, int n, int m)
{
    PyObject *pythonFloat, *pythonMatrix, *pythonList;
    int i, j;

    pythonMatrix = PyList_New(n);
    for(i = 0; i < n; i++)
    {
        pythonList = PyList_New(m);
        for(j = 0; j < m; j++)
        {
            pythonFloat = Py_BuildValue("f", matrix[i][j]);
            PyList_SetItem(pythonList, j, pythonFloat);
        }
        PyList_SetItem(pythonMatrix, i, pythonList);
    }
    freeMatrix(matrix, n);
    return pythonMatrix;
}


double** PyobjectToMatrix(PyObject *obj, int n)
{
    PyObject *currObj;
    double **matrix;
    int m;

    matrix = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        matrix[i] = (double*)calloc(n, sizeof(double));
        currObj = PyList_GetItem(obj, i);
        if (!PyList_Check(currObj))
        {
            PyErr_SetString(PyExc_TypeError, "The input must be a list of lists");
            return NULL;
        }
        m = (int)PyList_Size(currObj);
        for (int j = 0; j < m; j++)
        {
            matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(currObj, j));
        }
    }
    return matrix;
}


static PyObject* symnmf(PyObject *self, PyObject *args)
{
    PyObject *H;
    PyObject *W;


    if(!PyArg_ParseTuple(args, "OO", &H, &W)) {
        return NULL;
    }

    double **HMatrix = PyobjectToMatrix(H, (int)PyList_Size(H));
    double **WMatrix = PyobjectToMatrix(W, (int)PyList_Size(W));
    if (HMatrix == NULL || WMatrix == NULL) return NULL;

    int n = (int)PyList_Size(H);
    int m = (int)PyList_Size(PyList_GetItem(H, 1));
    double **next = initializeMatrix(n, m);
    double **HMatrixT = initializeMatrix(m, n);
    double **WHMatrix = initializeMatrix(n, m);
    double **HHtMatrix = initializeMatrix(n, n);
    double **HHtHMatrix = initializeMatrix(n, m);

    int i;
    int j;
    int iter = 0;
    double distance = INFINITY;
    double calc;

    while (iter < MAX_ITERATIONS)
    {
        HMatrixT = transpose(HMatrix);
        WHMatrix = matrixMultiply(WMatrix, HMatrix, n, m);
        HHtMatrix = matrixMultiply(HMatrix, HMatrixT, n, n);
        HHtHMatrix = matrixMultiply(HHtMatrix, HMatrix, n, m);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                calc = HMatrix[i][j] * (1 - BETA + BETA * (WHMatrix[i][j] / HHtHMatrix[i][j]));
                next[i][j] = calc;
            }
        }
        distance = SquaredFrobeniusNorm(next, HMatrix);
        if (distance < EPSILON) break;
        HMatrix = next;
        iter++;
    }
    
    freeMatrix(HMatrix, n);
    freeMatrix(WMatrix, n);
    freeMatrix(WHMatrix, n);
    freeMatrix(HHtMatrix, n);
    freeMatrix(HHtHMatrix, n);
    return matrixToPyObject(HMatrix, n, m);
}


static PyObject* sym(PyObject *self, PyObject *args)
{
    PyObject *X, retMat;
    struct vector* points;
    int numOfPoints;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL;
    }
    points = PyObjectToLinkedList(X);
    numOfPoints = countVectors(points);
    double** similarityMat = similarityMatrix(points, numOfPoints);
    retMat = matrixToPyObject(similarityMat, numOfPoints, numOfPoints);

    freeVector(points);
    freeMatrix(similarityMat, numOfPoints);    
    return retMat;
}


static PyObject* ddg(PyObject *self, PyObject *args)
{
    PyObject *X, *retMat;
    struct vector* points;
    int numOfPoints;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL;
    }
    points = PyObjectToLinkedList(X);
    numOfPoints = countVectors(points);
    double** similarityMat = similarityMatrix(points, numOfPoints);
    double** diagonalDegreeMat = diagonalDegreeMatrix(similarityMat, numOfPoints);
    retMat = matrixToPyObject(diagonalDegreeMat, numOfPoints, numOfPoints);

    freeVector(points);
    freeMatrix(similarityMat, numOfPoints);
    freeMatrix(diagonalDegreeMat, numOfPoints);
    return retMat;
}

static PyObject* norm(PyObject *self, PyObject *args)
{
    PyObject *X, *retMat;
    struct vector* points;
    int numOfPoints;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL;
    }
    points = PyObjectToLinkedList(X);
    numOfPoints = countVectors(points);
    double** similarityMat = similarityMatrix(points, numOfPoints);
    double** diagonalDegreeMat = diagonalDegreeMatrix(similarityMat, numOfPoints);
    double** normalizedSimilarityMat = normalizedSimilarityMatrix(similarityMat, diagonalDegreeMat, numOfPoints);
    retMat = matrixToPyObject(normalizedSimilarityMat, numOfPoints, numOfPoints);

    freeVector(points);
    freeMatrix(similarityMat, numOfPoints);
    freeMatrix(diagonalDegreeMat, numOfPoints);
    freeMatrix(normalizedSimilarityMat, numOfPoints);
    return retMat;
}

static PyMethodDef symnmfMethods[] = {
    {"symnmf",  
      (PyCFunction) symnmf,  
      METH_VARARGS,           
      PyDoc_STR("Perform the symnmf clustering algorithm")},
    {"sym",  
      (PyCFunction) sym,  
      METH_VARARGS,           
      PyDoc_STR("Calculate and output the similarity matrix")},
    {"ddg",  
      (PyCFunction) ddg,  
      METH_VARARGS,           
      PyDoc_STR("Calculate and output the Diagonal Degree Matrix")},
    {"norm",  
      (PyCFunction) norm,  
      METH_VARARGS,           
      PyDoc_STR("Calculate and output the normalized similarity matrix")},
    {NULL, NULL, 0, NULL} 
};


static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf", 
    NULL, 
    -1,
    symnmfMethods 
};


PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}
