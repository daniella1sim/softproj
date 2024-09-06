#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "symnmf.h"

#define EPSILON 0.0001
#define MAX_ITERATIONS 300
#define BETA 0.5


/**
 * @brief Print a linked list of coordinates
 * 
 * The function prints a linked list of coordinates.
 * 
 * @param c - A pointer to the head of the linked list 
 */
void printCord(struct cord *c)
{
    while(c != NULL)
    {
        printf("%.4f, ", c->value);
        c = c->next;
    }
}


/**
 * @brief Print a linked list of vectors
 * 
 * The function prints a linked list of vectors.
 * 
 * @param v - A pointer to the head of the linked list
 */
void printVector(struct vector *v)
{
    while(v != NULL)
    {
        printCord(v->cords);
        printf("\n");
        v = v->next;
    }
}


/**
 * @brief Convert a Python list of lists to a linked list of vectors
 * 
 * The function converts a Python list of lists (where each inner list represents a vector) 
 * to a linked list of `vector` structures.
 * 
 * @param obj - A Python object representing a list of lists.
 * @return A pointer to the head of the linked list of vectors, or NULL if an error occurs.
 */
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


/**
 * @brief Convert a matrix to a Python object
 * 
 * The function converts a `Matrix` structure to a Python list of lists.
 * 
 * @param matrix - A pointer to the `Matrix` structure to be converted.
 * @return A Python object representing the matrix, or NULL if an error occurs.
 */
PyObject* matrixToPyObject(Matrix *matrix)
{
    PyObject *pythonFloat;
    PyObject *pythonMatrix;
    PyObject *pythonList;
    int i;
    int j;

    pythonMatrix = PyList_New(matrix->rows);
    for(i = 0; i < matrix->rows; i++)
    {
        pythonList = PyList_New(matrix->cols);
        for(j = 0; j < matrix->cols; j++)
        {
            pythonFloat = Py_BuildValue("f", matrix->data[i][j]);
            PyList_SetItem(pythonList, j, pythonFloat);
        }
        PyList_SetItem(pythonMatrix, i, pythonList);
    }
    return pythonMatrix;
}


/**
 * @brief Convert a Python list of lists to a matrix
 * 
 * The function converts a Python list of lists (where each inner list represents a row of the matrix) 
 * to a `Matrix` structure.
 * 
 * @param obj - A Python object representing a list of lists.
 * @return A pointer to the `Matrix` structure, or NULL if an error occurs.
 */
Matrix *PyobjectToMatrix(PyObject *obj)
{
    PyObject *currObj;
    Matrix *matrix;
    int n;
    int m;
    int i;
    int j;
    
    if (!PyList_Check(obj))
    {
        PyErr_SetString(PyExc_TypeError, "The input must be a list of lists");
        return NULL;
    }
    n = (int)PyList_Size(obj);
    if (n == 0)
    {
        PyErr_SetString(PyExc_ValueError, "The input list is empty");
        return NULL;
    }
    
    currObj = PyList_GetItem(obj, 0);
    if (!PyList_Check(currObj)) {
        PyErr_SetString(PyExc_TypeError, "The input list must contain lists");
        return NULL;
    }
    m = (int)PyList_Size(currObj);
    
    matrix = initializeMatrix(n, m);
    
    for (i = 0; i < n; i++)
    {   
        currObj = PyList_GetItem(obj, i);
        if (!PyList_Check(currObj))
        {
            PyErr_SetString(PyExc_TypeError, "The input must be a list of lists");
            freeMatrix(matrix);
            return NULL;
        }
        for (j = 0; j < m; j++)
        {
            matrix->data[i][j] = PyFloat_AsDouble(PyList_GetItem(currObj, j));
        }
    }
    return matrix;
}


/**
 * @brief Perform Symmetric Non-negative Matrix Factorization (SymNMF)
 * 
 * The function performs the Symmetric Non-negative Matrix Factorization algorithm on 
 * the given matrices H and W. It iteratively updates the matrices to minimize the 
 * reconstruction error.
 * 
 * @param self Pointer to the module object (not used).
 * @param args A tuple containing the two matrices H and W as Python objects.
 * @return A Python object representing the factorized matrix, or NULL if an error occurs.
 */
static PyObject* symnmf(PyObject *self, PyObject *args)
{
    PyObject *H;
    PyObject *W;

    if(!PyArg_ParseTuple(args, "OO", &H, &W)) {
        return NULL;
    }

    Matrix *HMatrix = PyobjectToMatrix(H);
    Matrix *WMatrix = PyobjectToMatrix(W);
    if (HMatrix == NULL || WMatrix == NULL) return NULL;

    int n = HMatrix->rows;
    int m = HMatrix->cols;
    Matrix *next = initializeMatrix(n, m);
    Matrix *HMatrixT = initializeMatrix(m, n);
    Matrix *WHMatrix = initializeMatrix(n, m);
    Matrix *HHtMatrix = initializeMatrix(n, n);
    Matrix *HHtHMatrix = initializeMatrix(n, m);

    int i;
    int j;
    int iter = 0;
    double distance = INFINITY;
    double calc;
    
    while (iter < MAX_ITERATIONS)
    {
        Matrix* tempmatrix = transpose(HMatrix);
        freeMatrix(HMatrixT);
        HMatrixT = tempmatrix;
        
        tempmatrix = matrixMultiply(WMatrix, HMatrix);
        freeMatrix(WHMatrix);
        WHMatrix = tempmatrix;

        tempmatrix = matrixMultiply(HMatrix, HMatrixT);
        freeMatrix(HHtMatrix);
        HHtMatrix = tempmatrix;

        tempmatrix = matrixMultiply(HHtMatrix, HMatrix);
        freeMatrix(HHtHMatrix);
        HHtHMatrix = tempmatrix;

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                if (HHtHMatrix->data[i][j] != 0)
                {
                    calc = HMatrix->data[i][j] * (1 - BETA + BETA * (WHMatrix->data[i][j] / HHtHMatrix->data[i][j]));
                } else 
                {
                    calc = HMatrix->data[i][j] * (1 - BETA); // Handle division by zero
                }
                next->data[i][j] = calc;
            }
        }
        distance = MatrixDistance(next, HMatrix);
        if (distance < EPSILON) break;
        Matrix *tmpmatrixNext = next;
        next = HMatrix;
        HMatrix = tmpmatrixNext;
        

        iter++;
    }
    
    freeMatrix(HMatrix);
    freeMatrix(WMatrix);
    freeMatrix(WHMatrix);
    freeMatrix(HHtMatrix);
    freeMatrix(HHtHMatrix);
    
    PyObject *retMat = matrixToPyObject(next);
    freeMatrix(next);
    
    return retMat;
}


/**
 * @brief Calculate the similarity matrix
 * 
 * The function calculates the similarity matrix for a set of points.
 * 
 * @param self Pointer to the module object (not used).
 * @param args A tuple containing the points as a Python object.
 * @return A Python object representing the similarity matrix, or NULL if an error occurs.
 */
static PyObject* sym(PyObject *self, PyObject *args)
{
    PyObject *X, *retMat;
    struct vector* points;
    int numOfPoints;

    if(!PyArg_ParseTuple(args, "O", &X)) {
        return NULL;
    }
    points = PyObjectToLinkedList(X);
    numOfPoints = countVectors(points);
    Matrix *similarityMat = similarityMatrix(points, numOfPoints);
    retMat = matrixToPyObject(similarityMat);

    freeVector(points);
    freeMatrix(similarityMat);    
    return retMat;
}


/**
 * @brief Calculate the Diagonal Degree Matrix
 * 
 * The function calculates the Diagonal Degree Matrix from a set of points.
 * 
 * @param self Pointer to the module object (not used).
 * @param args A tuple containing the points as a Python object.
 * @return A Python object representing the Diagonal Degree Matrix, or NULL if an error occurs.
 */
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
    Matrix *similarityMat = similarityMatrix(points, numOfPoints);
    Matrix *diagonalDegreeMat = diagonalDegreeMatrix(similarityMat);
    retMat = matrixToPyObject(diagonalDegreeMat);
    freeVector(points);
    freeMatrix(similarityMat);
    freeMatrix(diagonalDegreeMat);
    return retMat;
}


/**
 * @brief Calculate the normalized similarity matrix
 * 
 * The function calculates the normalized similarity matrix from a set of points.
 * 
 * @param self Pointer to the module object (not used).
 * @param args A tuple containing the points as a Python object.
 * @return A Python object representing the normalized similarity matrix, or NULL if an error occurs.
 */
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
    Matrix *similarityMat = similarityMatrix(points, numOfPoints);
    Matrix *diagonalDegreeMat = diagonalDegreeMatrix(similarityMat);
    Matrix *normalizedSimilarityMat = normalizedSimilarityMatrix(similarityMat, diagonalDegreeMat);
    retMat = matrixToPyObject(normalizedSimilarityMat);

    freeVector(points);
    freeMatrix(similarityMat);
    freeMatrix(diagonalDegreeMat);
    freeMatrix(normalizedSimilarityMat);
    return retMat;
}


/* Module method definitions */
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


/* Module definition */
static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf", 
    NULL, 
    -1,
    symnmfMethods 
};


/**
 * @brief Initialize the symnmf module
 * 
 * This function is called when the module is imported. It creates the module object.
 * 
 * @return A pointer to the module object, or NULL if an error occurs.
 */
PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}
