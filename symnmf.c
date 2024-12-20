#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int i;
int j;
int k;


/**
 * @brief Struct for a single coordinate
 * 
 * A single coordinate is a value in a vector.
 */
struct cord{
    double value;
    struct cord *next;
};


/**
 * @brief Struct for a vector
 * 
 * A vector is a list of coordinates.
 */
struct vector{
    struct vector *next;
    struct cord *cords;
};


/**
 * @brief Struct for a 2D matrix
 * 
 * A 2D matrix is a list of lists of double values.
 */
typedef struct{
    double **data;
    int rows;
    int cols;
} Matrix;


/**
 * @brief Free dynamically allocated matrix
 * 
 * The function frees a dynamically allocated matrix.
 * 
 * @param matrix - A pointer to a Matrix struct
 */
void freeMatrix(Matrix *matrix){
    if (matrix == NULL) return;
    for (i = 0; i < matrix->rows; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);
    free(matrix);
}



/**
 * @brief Initialize a 2D matrix
 * 
 * The function allocates memory for a 2D array and returns a Matrix struct.
 * 
 * @param n - The number of rows in the matrix
 * @param m - The number of columns in the matrix
 * @return Matrix with all values initialized to 0
 */
Matrix *initializeMatrix(int n, int m)
{
    Matrix *mat = (Matrix*)calloc(1, sizeof(Matrix));
    if (mat == NULL) return NULL;
    mat->rows = n;
    mat->cols = m;
    mat->data = (double**)calloc(n, sizeof(double*));
    if (mat->data == NULL){
        freeMatrix(mat);
        return NULL;
    }
    for (i = 0; i < n; i++){
        mat->data[i] = (double*)calloc(m, sizeof(double));
        if (mat->data[i] == NULL){
            freeMatrix(mat);
            return NULL;
        };
    }
    return mat;
}


/**
 * @brief Transpose a given matrix
 * 
 * The function transposes a given matrix.
 * 
 * @param matrix - A pointer to a Matrix struct
 * @return Matrix - A transposed Matrix struct
 */
Matrix *transpose(Matrix *matrix){
    Matrix *transposed = initializeMatrix(matrix->cols, matrix->rows);
    if (transposed == NULL) return NULL;
    for (i = 0; i < matrix->cols; i++){
        for (j = 0; j < matrix->rows; j++){
            transposed->data[i][j] = matrix->data[j][i];
        }
    }
    return transposed;
}


/**
 * @brief Calculate the distance between two matrices
 * 
 * Calculate the distance between two matrices.
 * 
 * @param matrixA - The first Matrix struct
 * @param matrixB - The second Matrix struct
 * @return double - The distance between the two matrices
 */
double MatrixDistance(Matrix *matrixA, Matrix *matrixB){
    double sum = 0;
    if (matrixA->rows != matrixB->rows || matrixA->cols != matrixB->cols) {
        return -1;
    }

    for (i = 0; i < matrixA->rows; i++){
        for (j = 0; j < matrixA->cols; j++){  
            sum += (matrixA->data[i][j] - matrixB->data[i][j]) * (matrixA->data[i][j] - matrixB->data[i][j]);
        }
    }
    return sum;
}


/**
 * @brief Create a new coordinate
 * 
 * The function allocates memory for a new coordinate and returns a pointer to its head.
 * 
 * @return cord* 
 */
struct cord* createNewCord(){
    struct cord *headCord;
    headCord = (struct cord*)calloc(1, sizeof(struct cord));
    if (headCord == NULL) return NULL;
    headCord->next = NULL;
    return headCord;
}


/**
 * @brief Create a new vector
 * 
 * The function allocates memory for a new vector and returns a pointer to its head.
 * 
 * @return vector* 
 */
struct vector* createNewVector(){
    struct vector *headVector;
    headVector = (struct vector*)calloc(1, sizeof(struct vector));
    if (headVector == NULL) return NULL;
    headVector->next = NULL;
    return headVector;
}


/**
 * @brief Add a new coordinate
 * 
 * The function allocates memory for a new coordinate and returns a pointer to it.
 * 
 * @param currCord - A pointer to the current coordinate to add
 * @return cord* 
 */
struct cord* addNewCord(struct cord *currCord){
    currCord->next = (struct cord*)calloc(1, sizeof(struct cord));
    if (currCord->next == NULL) return NULL;
    currCord = currCord->next;
    currCord->next = NULL;
    return currCord;
}


/**
 * @brief Add a new vector
 * 
 * The function allocates memory for a new vector and returns a pointer to it.
 * 
 * @param currVec - A pointer to the current vector to add
 * @return struct vector* - A pointer to the new vector
 */
struct vector* addVector(struct vector *currVec){
    currVec->next = (struct vector*)calloc(1, sizeof(struct vector));
    if (currVec->next == NULL) return NULL;
    currVec = currVec->next;
    currVec->next = NULL;
    return currVec;
}


/**
 * @brief Free a linked list of coordinates
 * 
 * The function frees the memory allocated for a linked list of coordinates.
 * 
 * @param c - A pointer to the head of the linked list
 */
void freeCord(struct cord *c){
    struct cord *tmp;
    tmp = c;
    while (c != NULL){
        tmp = c;
        c = c->next;
        free(tmp);
    }
    free(c);
}


/**
 * @brief Free a linked list of vectors
 * 
 * The function frees the memory allocated for a linked list of vectors.
 * 
 * @param v - A pointer to the head of the linked list
 */
void freeVector(struct vector *v){
    struct vector *tmp;
    tmp = v;
    while (v->next != NULL){
        tmp = v;
        v = v->next;
        freeCord(tmp->cords);
        free(tmp);
    }
    freeCord(v->cords);
    free(v);
}


/**
 * @brief Load points from file
 * 
 * The function reads the file and loads the points into a linked list of vectors.
 * Each vector contains a linked list of cords.
 * 
 * @param file
 * @return struct vector* 
 */
struct vector* loadPoints(FILE *file){
    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    double n;
    char c;
    head_cord = (struct cord*)calloc(1, sizeof(struct cord));
    curr_cord = head_cord;
    curr_cord->next = NULL;
    head_vec = (struct vector*)calloc(1,sizeof(struct vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;

    while (fscanf(file, "%lf%c", &n, &c) == 2){
        if (c == '\n'){
            curr_cord->value = n;
            curr_vec->cords = head_cord;
            curr_vec->next = (struct vector*)calloc(1, sizeof(struct vector));
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = (struct cord*)calloc(1, sizeof(struct cord));
            curr_cord = head_cord;
            curr_cord->next = NULL;
            continue;
        }
        curr_cord->value = n;
        curr_cord->next = (struct cord*)calloc(1, sizeof(struct cord));
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
    }
    freeCord(head_cord);
    return head_vec;
}


/**
 * @brief Count the number of vectors
 * 
 * The function counts the number of vectors in the linked list.
 * 
 * @param headVec
 * @return int 
 */
int countVectors(struct vector *headVec){
    struct vector *currVec;
    int counter;

    currVec = headVec;
    counter = 0;
    while(currVec != NULL){
        counter++;
        currVec = currVec->next;
    }
    return counter - 1;
}


/**
 * @brief Compare two strings
 * 
 * The function compares two strings and returns 1 if they are equal and 0 otherwise.
 * 
 * @param s1 - The first string
 * @param s2 - The second string
 * @return int 
 */
int compareStrings(const char *s1, const char *s2){
    while(*s1 != '\0'){
        if (*s1 != *s2) return 0;
        s1++, s2++;
    }
    if (*s2 == '\0') return 1;
    return 0;
}


/**
 * @brief Calculate the distance between two points
 * 
 * The function calculates the Gaussian distance between two points.
 *
 * @param point1 - A pointer to the head of a linked list of cords representing the first point
 * @param point2 - A pointer to the head of a linked list of cords representing the second point
 * @return double - The distance between the two points
 */

double distance(struct cord *point1, struct cord *point2){
    double sum = 0;
    double val1, val2;

    while (point1 != NULL){
        val1 = point1->value;
        val2 = point2->value;
        sum += (val1 - val2) * (val1 - val2);
        point1 = point1->next;
        point2 = point2->next;
    }
    return exp(-sum / 2);
}


/**
 * @brief Calculate the sum of a column in a matrix
 * 
 * The function calculates the sum of a column in a matrix.
 * 
 * @param matrix - A pointer to a Matrix struct
 * @param col - The column to sum
 * @return double - The sum of the column
 */
double columnSum(Matrix *matrix, int col){
    double sum = 0;
    for (j = 0; j < matrix->rows; j++){
        sum += matrix->data[j][col];
    }
    return sum;
}


/**
 * @brief Multiply two matrices
 * 
 * The function multiplies two matrices.
 * 
 * @param matrix1 - The first Matrix struct
 * @param matrix2 - The second Matrix struct
 * @return Matrix - A Matrix struct representing the product of the two matrices
 */
Matrix *matrixMultiply(Matrix *matrix1, Matrix *matrix2){
    Matrix *result;
    if (matrix1->cols != matrix2->rows) {
        return NULL;
    }

    result = initializeMatrix(matrix1->rows, matrix2->cols);
    if (result == NULL) return NULL;

    for (i = 0; i < matrix1->rows; i++){
        for (j = 0; j < matrix2->cols; j++){
            result->data[i][j] = 0;
            for (k = 0; k < matrix1->cols; k++){
                result->data[i][j] += matrix1->data[i][k] * matrix2->data[k][j];
            }
        }
    }
    return result;
}


/**
 * @brief Print a matrix
 * 
 * The function prints a matrix.
 * 
 * @param matrix - A Matrix struct
 */
void printMatrix(Matrix *matrix){
    for(i = 0; i < matrix->rows; i++){
        for (j = 0; j < matrix->cols - 1; j++){
            printf("%.4f,", matrix->data[i][j]);
        }
        printf("%.4f\n", matrix->data[i][matrix->cols - 1]);
    }
}


/**
 * @brief Calculate the similarity matrix
 * 
 * The function calculates the similarity matrix for the given points.
 * 
 * @param points - A pointer to the head of a linked list of vector structures. Each vector represents a point.
 * @param numOfPoints - The number of points
 * @return Matrix - A Matrix struct representing the similarity matrix
 * 
 */
Matrix *similarityMatrix(struct vector *points, int numOfPoints){
    struct vector *currPointRow;
    struct vector *currPointColumn;
    int row;
    int column;
    Matrix *matrix;
    matrix = initializeMatrix(numOfPoints, numOfPoints);
    if (matrix == NULL) return NULL;
    currPointRow = points;
    for (row = 0; row < numOfPoints; row++){
        currPointColumn = points;
        for (column = 0; column < numOfPoints; column++){
            if (row == column) matrix->data[row][column] = 0;
            else{
                matrix->data[row][column] = distance(currPointColumn->cords, currPointRow->cords);
            }
            currPointColumn = currPointColumn->next;
        }       
        currPointRow = currPointRow->next;
    }
    return matrix;
}


/**
 * @brief Calculate the diagonal degree matrix
 * 
 * The function calculates the diagonal degree matrix for the given similarity matrix.
 * 
 * @param similarityMatrix - A pointer to the similarity matrix
 * @return Matrix - A Matrix struct representing the diagonal degree matrix
 */
Matrix *diagonalDegreeMatrix(Matrix *similarityMatrix){
    Matrix *ddg = NULL;
    double sum;
    ddg = initializeMatrix(similarityMatrix->rows, similarityMatrix->rows);
    if (ddg == NULL) return NULL;
    for (i = 0; i < similarityMatrix->rows; i++){
        sum = columnSum(similarityMatrix, i);
        ddg->data[i][i] = sum;
    }
    return ddg;
}


/**
 * @brief Calculate the inverse square root of the diagonal degree matrix
 * 
 * The function calculates the inverse square root of the diagonal degree matrix.
 * 
 * @param diagonalDegreeMatrix - A pointer to the diagonal degree matrix
 * @return Matrix - A Matrix struct representing the inverse square root of the diagonal degree matrix
 */
Matrix *diagInvSqrtMatrix(Matrix *diagonalDegreeMatrix){
    Matrix *norm = initializeMatrix(diagonalDegreeMatrix->rows, diagonalDegreeMatrix->cols);
    if (norm == NULL) return NULL;
    for (i = 0; i < diagonalDegreeMatrix->rows; i++){
        if (diagonalDegreeMatrix->data[i][i] != 0) {
            norm->data[i][i] = 1 / sqrt(diagonalDegreeMatrix->data[i][i]);
        }
        else norm->data[i][i] = 0;
    }
    return norm;
}


/**
 * @brief Calculate the normalized similarity matrix
 * 
 * The function calculates the normalized similarity matrix for the given similarity matrix and diagonal degree matrix.
 * 
 * @param similarityMat - A pointer to the simmilarity matrix
 * @param diagonalDegreeMat - A pointer to the diagonal degree matrix
 * @return Matrix - A Matrix struct representing the normalized similarity matrix
 */
Matrix *normalizedSimilarityMatrix(Matrix *similarityMat, Matrix *diagonalDegreeMat){
    Matrix *diagInvSqrtMat, *first, *result;
    diagInvSqrtMat = diagInvSqrtMatrix(diagonalDegreeMat);
    if (diagInvSqrtMat == NULL) return NULL;
    first = matrixMultiply(diagInvSqrtMat, similarityMat);
    if (first == NULL){
        freeMatrix(diagInvSqrtMat);
        return NULL;
    };
    result = matrixMultiply(first, diagInvSqrtMat);
    if (result == NULL){
        freeMatrix(diagInvSqrtMat);
        freeMatrix(first);
        return NULL;
    }    
    freeMatrix(diagInvSqrtMat);
    freeMatrix(first);
    return result;
}


/**
 * @brief Load points from a file and return a linked list of vectors.
 * 
 * This function opens a file, reads the points using the loadPoints function, and then 
 * closes the file. It returns the head of the linked list of points (vectors).
 * 
 * @param fileName Name of the file to read the points from.
 * @return struct vector* Pointer to the head of the linked list of vectors, or NULL on failure.
 */
struct vector* loadPointsFromFile(const char *fileName) {
    struct vector *points;
    FILE *file = fopen(fileName, "r");
    if (!file) return NULL;

    points = loadPoints(file);  /* Load points from the file */
    fclose(file);  /* Close the file */
    return points;
}


/**
 * @brief Process the goal provided by the user and compute the corresponding matrix.
 * 
 * This function checks the user's goal (sym, ddg, norm) and computes the corresponding matrix
 * (similarity matrix, diagonal degree matrix, or normalized similarity matrix). It also prints
 * the matrix based on the goal and frees any intermediate matrices.
 * 
 * @param goal The goal provided by the user (sym, ddg, norm).
 * @param similarityMat The similarity matrix computed from the points.
 */
int processGoal(const char *goal, Matrix *similarityMat) {
    char sym[] = "sym";
    char ddg[] = "ddg";
    char norm[] = "norm";
    Matrix *diagonalDegreeMat = NULL;
    Matrix *normalizedSimilarityMat = NULL;

    if (compareStrings(goal, sym) == 1) printMatrix(similarityMat);

    else {
        diagonalDegreeMat = diagonalDegreeMatrix(similarityMat);
        if (diagonalDegreeMat == NULL) return -1;
        if (compareStrings(goal, ddg) == 1) {
            printMatrix(diagonalDegreeMat);
            freeMatrix(diagonalDegreeMat);
        } 
        else if (compareStrings(goal, norm) == 1) {
            normalizedSimilarityMat = normalizedSimilarityMatrix(similarityMat, diagonalDegreeMat);
            if (normalizedSimilarityMat == NULL) {
                freeMatrix(diagonalDegreeMat);
                return -1;
            }
            printMatrix(normalizedSimilarityMat);
            freeMatrix(diagonalDegreeMat);
            freeMatrix(normalizedSimilarityMat);
        } 
        else printf("Invalid goal: %s\n", goal);
    } 
    return 0;
}


/**
 * @brief Main function
 * 
 * The main function reads the goal and the file name from the command line arguments.
 * It then calculates the similarity matrix, the diagonal degree matrix, or the normalized similarity matrix based on the goal.
 * After the calculation, the function prints the matrix and frees the memory.
 * 
 * @param argc - The number of command line arguments
 * 
 */
int main(int argc, char *argv[]) {
    char *goal, *fileName;
    struct vector *points;
    int numOfPoints, err;
    Matrix *similarityMat;

    if (argc != 3){
        printf("An Error Has Occured!\n");
        return 1;
    }
    
    goal = argv[1];
    fileName = argv[2];
    points = loadPointsFromFile(fileName);
    if (points == NULL){
        printf("An Error Has Occured!\n");
        return 1;
    }

    numOfPoints = countVectors(points);
    similarityMat = similarityMatrix(points, numOfPoints);
    if (similarityMat == NULL) {
        freeVector(points);
        return 1;
    }
    err = processGoal(goal, similarityMat);
    freeMatrix(similarityMat);
    freeVector(points);
    return (err == -1) ? 1 : 0;
}
