# define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "symnmf.h"

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

PyObject* matrixToPyObject(double** matrix, int n)
{
    PyObject *pythonFloat, *pythonMatrix, *pythonList;
    int i, j;

    pythonMatrix = PyList_New(n);
    for(i = 0; i < n; i++)
    {
        pythonList = PyList_New(n);
        for(j = 0; j < n; j++)
        {
            pythonFloat = Py_BuildValue("f", matrix[i][j]);
            PyList_SetItem(pythonList, j, pythonFloat);
        }
        PyList_SetItem(pythonMatrix, i, pythonList);
    }
    return pythonMatrix;
}

static PyObject* symnmf(PyObject *self, PyObject *args)
{
    PyObject *H;
    PyObject *W;

    if(!PyArg_ParseTuple(args, "OO", &H, &W)) {
        return NULL;
    }
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
    retMat = matrixToPyObject(similarityMat, numOfPoints);

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
    retMat = matrixToPyObject(diagonalDegreeMat, numOfPoints);

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
    retMat = matrixToPyObject(normalizedSimilarityMat, numOfPoints);

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