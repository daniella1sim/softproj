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
 * @brief Get the size of the specified row in the list of lists.
 * 
 * This function retrieves the number of elements in a specified inner list
 * to determine the number of columns for the vectors.
 * 
 * @param obj A Python object representing a list of lists.
 * @param rowIndex The index of the row to check.
 * @return The size of the row, or -1 if an error occurs.
 */
int getColumnSize(PyObject *obj, int rowIndex) {
    PyObject *rowObj = PyList_GetItem(obj, rowIndex);
    if (!PyList_Check(rowObj)) {
        PyErr_SetString(PyExc_TypeError, "The input list must contain lists");
        return -1; 
    }
    return (int)PyList_Size(rowObj);
}


/**
 * @brief Populate a vector structure from a Python list.
 * 
 * This function populates the given vector using values from the corresponding 
 * inner list from the Python object.
 * 
 * @param vec Pointer to the vector structure being populated.
 * @param obj A Python object representing the inner list.
 * @param columns The number of elements in the vector.
 * @return int 1 if successful, 0 if an error occurs.
 */
int populateVector(struct vector *vec, PyObject *obj, int columns) {
    struct cord *headCord = createNewCord();
    struct cord *currCord = headCord;

    for (int j = 0; j < columns; j++) {
        currCord->value = PyFloat_AsDouble(PyList_GetItem(obj, j));
        if (j < columns - 1) {
            currCord = addNewCord(currCord);
        }
    }

    vec->cords = headCord; 
    return 1; 
}


/**
 * @brief Convert a Python list of lists to a linked list of vectors.
 * 
 * This function converts a Python list of lists (where each inner list represents a vector) 
 * to a linked list of `vector` structures.
 * 
 * @param obj A Python object representing a list of lists.
 * @return A pointer to the head of the linked list of vectors, or NULL if an error occurs.
 */
struct vector* PyObjectToLinkedList(PyObject *obj) {
    if (!PyList_Check(obj)) {
        PyErr_SetString(PyExc_TypeError, "Input must be a list of lists");
        return NULL;
    }

    int rows = (int)PyList_Size(obj);
    if (rows == 0) {
        return NULL;
    }

    int columns = getColumnSize(obj, 0);
    struct vector *headVec = createNewVector();
    struct vector *currVec = headVec;

    for (int i = 0; i < rows; i++) {
        PyObject *currObj = PyList_GetItem(obj, i);
        if (!populateVector(currVec, currObj, columns)) {
            freeLinkedList(headVec);
            return NULL;
        }
        if (i < rows - 1) {
            currVec = addVector(currVec);
        }
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
 * @brief Convert a Python list of lists to a matrix.
 * 
 * This function converts a Python list of lists (where each inner list represents a row of the matrix) 
 * to a `Matrix` structure.
 * 
 * @param obj A Python object representing a list of lists.
 * @return A pointer to the `Matrix` structure, or NULL if an error occurs.
 */
Matrix *PyobjectToMatrix(PyObject *obj) {
    if (!PyList_Check(obj)) {
        PyErr_SetString(PyExc_TypeError, "The input must be a list of lists");
        return NULL;
    }

    int n = (int)PyList_Size(obj);
    if (n == 0) {
        PyErr_SetString(PyExc_ValueError, "The input list is empty");
        return NULL;
    }

    int m;
    Matrix *matrix = initializeMatrix(n, getRowSize(obj, 0));
    if (matrix == NULL) {
        return NULL; // Memory allocation failed
    }

    for (int i = 0; i < n; i++) {
        if (!populateRow(matrix, obj, i, m)) {
            freeMatrix(matrix);
            return NULL; // Error occurred in populating row
        }
    }

    return matrix;
}



/**
 * @brief Update the intermediate matrices during SymNMF.
 * 
 * This function performs the necessary matrix transpositions and multiplications required
 * for updating the factor matrices during the SymNMF algorithm.
 * 
 * @param HMatrix Pointer to the H matrix.
 * @param WMatrix Pointer to the W matrix.
 * @param HMatrixT Pointer to the transposed H matrix.
 * @param WHMatrix Pointer to the matrix resulting from W * H.
 * @param HHtMatrix Pointer to the matrix resulting from H * H^T.
 * @param HHtHMatrix Pointer to the matrix resulting from H * H^T * H.
 */
void updateMatrices(Matrix *HMatrix, Matrix *WMatrix, Matrix **HMatrixT, Matrix **WHMatrix, 
                    Matrix **HHtMatrix, Matrix **HHtHMatrix) {
    Matrix* tempmatrix;

    /* Transpose HMatrix */
    tempmatrix = transpose(HMatrix);
    freeMatrix(*HMatrixT);
    *HMatrixT = tempmatrix;

    /* Multiply WMatrix and HMatrix */
    tempmatrix = matrixMultiply(WMatrix, HMatrix);
    *WHMatrix = tempmatrix;

    /* Multiply HMatrix and HMatrixT */
    tempmatrix = matrixMultiply(HMatrix, *HMatrixT);
    freeMatrix(*HHtMatrix);
    *HHtMatrix = tempmatrix;

    /* Multiply HHtMatrix and HMatrix */
    tempmatrix = matrixMultiply(*HHtMatrix, HMatrix);
    freeMatrix(*HHtHMatrix);
    *HHtHMatrix = tempmatrix;
}


/**
 * @brief Update the H matrix using WHMatrix and HHtHMatrix.
 * 
 * This function updates the elements of the H matrix using the calculated matrices
 * WHMatrix (W * H) and HHtHMatrix (H * H^T * H). It handles cases where division by zero
 * might occur.
 * 
 * @param HMatrix Pointer to the H matrix.
 * @param WHMatrix Pointer to the matrix W * H.
 * @param HHtHMatrix Pointer to the matrix H * H^T * H.
 * @param next Pointer to the matrix where the updated H values are stored.
 * @param n Number of rows in the matrices.
 * @param m Number of columns in the matrices.
 */
void updateHMatrix(Matrix *HMatrix, Matrix *WHMatrix, Matrix *HHtHMatrix, Matrix *next, int n, int m) {
    double calc;
    int i;
    int j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (HHtHMatrix->data[i][j] != 0) {
                calc = HMatrix->data[i][j] * (1 - BETA + BETA * (WHMatrix->data[i][j] / HHtHMatrix->data[i][j]));
            } else {
                calc = HMatrix->data[i][j] * (1 - BETA);  // Handle division by zero
            }
            next->data[i][j] = calc;
        }
    }
}



/**
 * @brief Perform SymNMF algorithm on matrices H and W.
 * 
 * This function iteratively updates the matrices H and W to minimize the reconstruction error.
 * It performs matrix transpositions and multiplications, updates the factor matrices, and checks 
 * for convergence using the Frobenius norm distance.
 * 
 * @param HMatrix Pointer to the H matrix.
 * @param WMatrix Pointer to the W matrix.
 * @return Matrix* Pointer to the resulting factorized matrix after convergence.
 */
Matrix* performSymNMF(Matrix *HMatrix, Matrix *WMatrix) {
    int n = HMatrix->rows;
    int m = HMatrix->cols;
    Matrix *next = initializeMatrix(n, m);
    Matrix *HMatrixT = initializeMatrix(m, n);
    Matrix *WHMatrix = initializeMatrix(n, m);
    Matrix *HHtMatrix = initializeMatrix(n, n);
    Matrix *HHtHMatrix = initializeMatrix(n, m);
    Matrix *tmpmatrixNext = next;
    int iter = 0;
    double distance = INFINITY;

    while (iter < MAX_ITERATIONS) {
        updateMatrices(HMatrix, WMatrix, &HMatrixT, &WHMatrix, &HHtMatrix, &HHtHMatrix);
        updateHMatrix(HMatrix, WHMatrix, HHtHMatrix, next, n, m);
        distance = MatrixDistance(next, HMatrix);
        if (distance < EPSILON) break;
        tmpmatrixNext = next;
        next = HMatrix;
        HMatrix = tmpmatrixNext;
        iter++;
    }

    freeMatrix(HMatrixT);
    freeMatrix(WHMatrix);
    freeMatrix(HHtMatrix);
    freeMatrix(HHtHMatrix);
    freeMatrix(tmpmatrixNext);
    return next;
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
    Matrix *HMatrix;
    Matrix *WMatrix;
    Matrix *result;
    PyObject *retMat;

    if(!PyArg_ParseTuple(args, "OO", &H, &W)) {
        return NULL;
    }

    HMatrix = PyobjectToMatrix(H);
    WMatrix = PyobjectToMatrix(W);
    if (HMatrix == NULL || WMatrix == NULL) return NULL;

    result = performSymNMF(HMatrix, WMatrix);
    retMat = matrixToPyObject(result);

    freeMatrix(result);
    
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
