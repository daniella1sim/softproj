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
 * @brief Frees 5 matrices
 * 
 */
void freeMatrices(Matrix *M1, Matrix *M2, Matrix *M3, Matrix *M4, Matrix *M5){
    freeMatrix(M1);
    freeMatrix(M2);
    freeMatrix(M3);
    freeMatrix(M4);
    freeMatrix(M5);
    return;
}


/**
 * @brief Print a linked list of coordinates
 * 
 * The function prints a linked list of coordinates.
 * 
 * @param c - A pointer to the head of the linked list 
 */
void printCord(struct cord *c){
    while(c != NULL){
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
void printVector(struct vector *v){
    while(v != NULL){
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
    struct cord *headCord; 
    struct cord *currCord;

    headCord = createNewCord();
    if (headCord == NULL) return -1;
    currCord = headCord;

    for (int j = 0; j < columns; j++) {
        currCord->value = PyFloat_AsDouble(PyList_GetItem(obj, j));
        if (j < columns - 1) {
            currCord = addNewCord(currCord);
            if (currCord == NULL) {
                freeCord(headCord);
                return -1;
            }
        }
    }

    vec->cords = headCord; 
    return 0; 
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
    struct vector *headVec;
    struct vector *currVec;
    int columns, rows, i, err;
    PyObject *currObj;

    if (!PyList_Check(obj)) return NULL;
    
    rows = (int)PyList_Size(obj);
    if (rows == 0) return NULL;

    columns = getColumnSize(obj, 0);
    headVec = createNewVector();
    if (headVec == NULL) return NULL;
    currVec = headVec;

    for (i = 0; i < rows; i++) {
        currObj = PyList_GetItem(obj, i);
        err = populateVector(currVec, currObj, columns);
        if (err == -1) return NULL;
        if (i < rows - 1) {
            currVec = addVector(currVec);
            if (currVec == NULL){
                freeVector(headVec);
                return NULL;
            }
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
 * @brief Populate a row of the matrix from a Python list.
 * 
 * This function populates a row of the matrix using values from the corresponding
 * 
 * inner list from the Python object.
 * 
 * @param matrix Pointer to the matrix structure being populated.
 * @param obj A Python object representing the list of lists.
 * @param row The index of the row being populated.
 * @param columns The number of columns in the matrix.
 * @return int 1 if successful, 0 if an error occurs.
 */
int populateRow(Matrix *matrix, PyObject *obj, int row, int columns) {
    int j;
    PyObject *rowObj;
    rowObj = PyList_GetItem(obj, row);
    
    if (!PyList_Check(rowObj)) {
        PyErr_SetString(PyExc_TypeError, "The input list must contain lists");
        return 0; 
    }
    for (j = 0; j < columns; j++) {
        matrix->data[row][j] = PyFloat_AsDouble(PyList_GetItem(rowObj, j));
    }
    return 1; 
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
    Matrix *matrix;
    int m, n, i;
    
    if (!PyList_Check(obj)) {
        PyErr_SetString(PyExc_TypeError, "The input must be a list of lists");
        return NULL;
    }

    n = (int)PyList_Size(obj);
    if (n == 0) {
        PyErr_SetString(PyExc_ValueError, "The input list is empty");
        return NULL;
    }
    m = getColumnSize(obj, 0);
    matrix = initializeMatrix(n, m);
    if (matrix == NULL) {
        return NULL; /* Memory allocation failed */
    }

    for (i = 0; i < n; i++) {
        if (!populateRow(matrix, obj, i, m)) {
            freeMatrix(matrix);
            return NULL; /* Error occurred in populating row */
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
int updateMatrices(Matrix *HMatrix, Matrix *WMatrix, Matrix **HMatrixT, Matrix **WHMatrix, Matrix **HHtMatrix, Matrix **HHtHMatrix) {
    Matrix* tempmatrix;
    tempmatrix = transpose(HMatrix);
    if (tempmatrix == NULL) return -1;
    freeMatrix(*HMatrixT);
    *HMatrixT = tempmatrix;
    tempmatrix = matrixMultiply(WMatrix, HMatrix);
    if (tempmatrix == NULL){
        freeMatrix(*HMatrixT);
        return -1;
    }
    *WHMatrix = tempmatrix;
    tempmatrix = matrixMultiply(HMatrix, *HMatrixT);
    if (tempmatrix == NULL){
        freeMatrices(*WHMatrix, *HMatrixT, NULL, NULL, NULL);
        return -1;
    }
    freeMatrix(*HHtMatrix);
    *HHtMatrix = tempmatrix;
    tempmatrix = matrixMultiply(*HHtMatrix, HMatrix);
    if (tempmatrix == NULL){
        freeMatrices(*WHMatrix, *HMatrixT, *HHtMatrix, NULL, NULL);
        return -1;
    } freeMatrix(*HHtHMatrix);
    *HHtHMatrix = tempmatrix;
    return 0;
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
                calc = HMatrix->data[i][j] * (1 - BETA);
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
    int err, iter;
    int n = HMatrix->rows;
    int m = HMatrix->cols;
    Matrix *next = initializeMatrix(n, m);
    Matrix *HMatrixT = initializeMatrix(m, n);
    Matrix *WHMatrix = initializeMatrix(n, m);
    Matrix *HHtMatrix = initializeMatrix(n, n);
    Matrix *HHtHMatrix = initializeMatrix(n, m);
    Matrix *tmpmatrixNext;
    if (next == NULL || HMatrixT == NULL || WHMatrix == NULL || HHtMatrix == NULL || HHtHMatrix == NULL) {
        freeMatrices(HMatrixT, WHMatrix, HHtMatrix, HHtHMatrix, next);
        return NULL;
    } tmpmatrixNext = next;
    iter = 0;
    double distance = INFINITY;
    while (iter < MAX_ITERATIONS) {
        err = updateMatrices(HMatrix, WMatrix, &HMatrixT, &WHMatrix, &HHtMatrix, &HHtHMatrix);
        if (err == -1) return NULL;
        updateHMatrix(HMatrix, WHMatrix, HHtHMatrix, next, n, m);
        distance = MatrixDistance(next, HMatrix);
        if (distance < EPSILON) break;
        tmpmatrixNext = next;
        next = HMatrix;
        HMatrix = tmpmatrixNext;
        iter++;
    } freeMatrices(HMatrixT, WHMatrix, HHtMatrix, HHtHMatrix, tmpmatrixNext);
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

    if(!PyArg_ParseTuple(args, "OO", &H, &W)) return NULL;

    HMatrix = PyobjectToMatrix(H);
    if (HMatrix == NULL) return NULL;
    WMatrix = PyobjectToMatrix(W);
    if (WMatrix == NULL) {
        freeMatrix(HMatrix);
        return NULL;
    }

    result = performSymNMF(HMatrix, WMatrix);
    if (result == NULL) {
        freeMatrices(result, HMatrix, WMatrix, NULL, NULL);
        return NULL;
    }
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
    if (points == NULL) return NULL;
    numOfPoints = countVectors(points) + 1;
    Matrix *similarityMat = similarityMatrix(points, numOfPoints);
    if (similarityMat == NULL) {
        freeVector(points);
        return NULL;
    }
    retMat = matrixToPyObject(similarityMat);
    if (retMat == NULL) {
        freeVector(points);
        freeMatrix(similarityMat);
        return NULL;
    }

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
    if (points == NULL) return NULL;
    numOfPoints = countVectors(points) + 1;
    Matrix *similarityMat = similarityMatrix(points, numOfPoints);
    if (similarityMat == NULL) {
        freeVector(points);
        return NULL;
    }
    Matrix *diagonalDegreeMat = diagonalDegreeMatrix(similarityMat);
    if (diagonalDegreeMat == NULL) {
        freeVector(points);
        freeMatrix(similarityMat);
        return NULL;
    }
    retMat = matrixToPyObject(diagonalDegreeMat);
    freeVector(points);
    freeMatrices(similarityMat, diagonalDegreeMat, NULL, NULL, NULL);
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
static PyObject* norm(PyObject *self, PyObject *args){
    PyObject *X, *retMat;
    struct vector* points;
    int numOfPoints;
    Matrix *similarityMat, *diagonalDegreeMat, *normalizedSimilarityMat;
    if(!PyArg_ParseTuple(args, "O", &X)) return NULL;
    points = PyObjectToLinkedList(X);
    if (points == NULL) return NULL;
    numOfPoints = countVectors(points) + 1;
    similarityMat = similarityMatrix(points, numOfPoints);
    if (similarityMat == NULL) {
        freeVector(points);
        return NULL;
    }
    diagonalDegreeMat = diagonalDegreeMatrix(similarityMat);
    if (diagonalDegreeMat == NULL) {
        freeVector(points);
        freeMatrix(similarityMat);
        return NULL;
    }
    normalizedSimilarityMat = normalizedSimilarityMatrix(similarityMat, diagonalDegreeMat);
    if (normalizedSimilarityMat == NULL) {
        freeVector(points);
        freeMatrices(similarityMat, diagonalDegreeMat, normalizedSimilarityMat, NULL, NULL);
        return NULL;
    }
    retMat = matrixToPyObject(normalizedSimilarityMat);
    freeVector(points);
    freeMatrices(similarityMat, diagonalDegreeMat, normalizedSimilarityMat, NULL, NULL);
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
