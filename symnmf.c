#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/**
 * @brief Struct for a single coordinate
 * 
 * A single coordinate is a value in a vector.
 */
struct cord
{
    double value;
    struct cord *next;
};


/**
 * @brief Struct for a vector
 * 
 * A vector is a list of coordinates.
 */
struct vector
{
    struct vector *next;
    struct cord *cords;
};

/**
 * @brief Initialize a 2D matrix
 * 
 * The function allocates memory for a 2D array and returns a pointer to its head.
 * 
 * @param n - The number of rows in the matrix
 * @param m - The number of columns in the matrix
 * @return double** 
 */
double ** initializeMatrix(int n, int m)
{
    double **matrix;
    matrix = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        matrix[i] = (double*)calloc(m, sizeof(double));
    }
    return matrix;
}


/**
 * @brief Transpose a given matrix
 * 
 * The function transposes a given matrix.
 * 
 * @param matrix - A pointer to the head of a 2D array
 * @return double** - A pointer to the head of a 2D array
 *  
 */
double **transpose(double **matrix)
{
    int n = sizeof(matrix) / sizeof(matrix[0]);
    int m = sizeof(matrix[0]) / sizeof(matrix[0][0]);
    double **transposed = initializeMatrix(m, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            transposed[j][i] = matrix[i][j];
        }
    }
    return transposed;
}


double SquaredFrobeniusNorm(double **matrixA, double** matrixB)
{
    double sum = 0;
    int n = sizeof(matrixA) / sizeof(matrixA[0]);
    if (n != sizeof(matrixB) / sizeof(matrixB[0])) return -1;
    int m = sizeof(matrixA[0]) / sizeof(matrixA[0][0]);
    if (m != sizeof(matrixB[0]) / sizeof(matrixB[0][0])) return -1;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {  
            sum += (matrixA[i][j] - matrixB[i][j]) * (matrixA[i][j] - matrixB[i][j]);
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
struct cord* createNewCord()
{
    struct cord *headCord;
    headCord = (struct cord*)calloc(1, sizeof(struct cord));
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
struct vector* createNewVector()
{
    struct vector *headVector;
    headVector = (struct vector*)calloc(1, sizeof(struct vector));
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
struct cord* addNewCord(struct cord *currCord)
{
    currCord->next = (struct cord*)calloc(1, sizeof(struct cord));
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
struct vector* addVector(struct vector *currVec)
{
    currVec->next = (struct vector*)calloc(1, sizeof(struct vector));
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
void freeCord(struct cord *c)
{
    struct cord *tmp;
    tmp = c;
    while (c->next != NULL)
    {
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
void freeVector(struct vector *v)
{
    struct vector *tmp;
    tmp = v;
    while (v->next != NULL)
    {
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
struct vector* loadPoints(FILE *file)
{
    char ch; // ch is used to check if the current character is a newline
    char prevCh; // prevCh is used to check if the previous character was a newline
    double currVal;

    struct vector *currVec; 
    struct vector *headVec;
    struct cord *headCord; 
    struct cord *currCord;

    headCord = (struct cord*)calloc(1, sizeof(struct cord));
    currCord = headCord;
    currCord->next = NULL;

    headVec = (struct vector*)calloc(1, sizeof(struct vector));
    currVec = headVec;
    currVec->next = NULL;

    prevCh = '\0';
    while (fscanf(file, "%lf%c", &currVal, &ch) == 2) // Read the value and the character
    {
        if(ch == '\n')
        {
            currCord->value = currVal;
            currVec->cords = headCord;
        }
        else
        {
            if (prevCh == '\n')
            {
                currVec = addVector(currVec);
                headCord = createNewCord();
                currCord = headCord;
            }
            currCord->value = currVal;
            currCord = addNewCord(currCord);
        }
        prevCh = ch;
    }
    return headVec;
}


/**
 * @brief Count the number of vectors
 * 
 * The function counts the number of vectors in the linked list.
 * 
 * @param headVec
 * @return int 
 */
int countVectors(struct vector *headVec)
{
    struct vector *currVec;
    int counter;

    currVec = headVec;
    counter = 0;
    while(currVec != NULL)
    {
        counter++;
        currVec = currVec->next;
    }
    return counter;
}


/**
 * @brief Free a 2D array
 * 
 * The function frees the memory allocated for a 2D array.
 * 
 * @param matrix - A pointer to the head of the 2D array
 * @param n - The number of rows in the matrix
 */
void freeMatrix(double **matrix, int n)
{
    for(int i = 0; i < n; i ++)
    {
        free(matrix[i]);
    }
    free(matrix);
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
int compareStrings(char *s1, char *s2)
{
    while(*s1 != '\0')
    {
        if (*s1 != *s2)
        {
            return 0;
        }
        s1++, s2++;
    }
    if (*s2 == '\0')
    {
        return 1;
    }
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

double distance(struct cord *point1, struct cord *point2)
{
    double sum;
    double val1, val2;

    sum = 0;
    while (point1 != NULL)
    {
        val1 = point1->value;
        val2 = point2->value;
        sum += ((val1 - val2)*(val1 - val2));
        point1 = point1->next;
        point2 = point2->next;
    }
    return exp(-sum/2);
}


/**
 * @brief Calculate the sum of a column in a matrix
 * 
 * The function calculates the sum of a column in a matrix.
 * 
 * @param matrix - A pointer to the head of a 2D array
 * @param col - The column to sum
 * @param n - The number of rows and columns in the matrix
 * @return double - The sum of the column
 */
double columnSum(double** matrix, int col, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += matrix[i][col];
    }
    return sum;
}


/**
 * @brief Multiply two matrices
 * 
 * The function multiplies two matrices.
 * 
 * @param matrix1 - A pointer to the head of the first 2D array
 * @param matrix2 - A pointer to the head of the second 2D array
 * @param n - The number of rows and columns in the matrices
 * @return double** - A pointer to the head of a 2D array representing the product of the two matrices
 */
double** matrixMultiply(double** matrix1, double** matrix2, int n, int m)
{
    double ij;
    double **matrix;
    matrix = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        matrix[i] = (double*)calloc(n, sizeof(double));
        for (int j = 0; j < n; j++)
        {
            ij = 0;
            for (int k = 0; k < m; k++)
            {
                ij += matrix1[i][k] * matrix2[k][j];
            }
            matrix[i][j] = ij;
        }
    }
    return matrix;
}

/**
 * @brief Print a matrix
 * 
 * The function prints a 2D array.
 * 
 * @param matrix - A pointer to the head of a 2D array
 * @param n - The number of rows and columns in the matrix
 */
void printMatrix(double** matrix, int n)
{
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            printf("%.4f,", matrix[i][j]);
        }
        printf("%.4f\n", matrix[i][n - 1]);
    }
}


/**
 * @brief Calculate the similarity matrix
 * 
 * The function calculates the similarity matrix for the given points.
 * 
 * @param points - A pointer to the head of a linked list of vector structures. Each vector represents a point.
 * @param numOfPoints - The number of points
 * @return double** - A pointer to the head of a 2D array representing the similarity matrix
 * 
 */
double** similarityMatrix(struct vector *points, int numOfPoints)
{
    struct vector *currPointRow;
    struct vector *currPointColumn;
    double **matrix; // The similarity matrix

    currPointRow = points;
    matrix = (double**)calloc(numOfPoints, sizeof(double*));
    for (int row = 0; row < numOfPoints; row++)
    {
        matrix[row] = (double*)calloc(numOfPoints, sizeof(double));
        currPointColumn = points;
        for (int column = 0; column < numOfPoints; column++)
        {
            if (row == column)
            {
                matrix[row][column] = 0;
            } else{
                matrix[row][column] = distance(currPointColumn->cords, currPointRow->cords);
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
 * @param similarityMatrix - A pointer to the head of a 2D array representing the similarity matrix
 * @param n - The number of rows and columns in the matrix
 * @return double** - A pointer to the head of a 2D array representing the diagonal degree matrix
 */
double** diagonalDegreeMatrix(double **similarityMatrix, int n)
{
    double **ddg;
    ddg = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        ddg[i] = (double*)calloc(n, sizeof(double));
        for (int j = 0; j < n; j++)
        {
            if (i == j) ddg[i][j] = columnSum(similarityMatrix, j, n);
            else ddg[i][j] = 0;
        }
    }
    return ddg;
}


/**
 * @brief Calculate the inverse square root of the diagonal degree matrix
 * 
 * The function calculates the inverse square root of the diagonal degree matrix.
 * 
 * @param diagonalDegreeMatrix - A pointer to the head of a 2D array representing the diagonal degree matrix
 * @param n - The number of rows and columns in the matrix
 * @return double** - A pointer to the head of a 2D array representing the inverse square root of the diagonal degree matrix
 */
double** diagInvSqrtMatrix(double **diagonalDegreeMatrix, int n)
{
    double **norm;
    norm = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        norm[i] = (double*)calloc(n, sizeof(double));
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                norm[i][j] = 1 / sqrt(diagonalDegreeMatrix[i][j]);
            } else {
                norm[i][j] = 0;
            }
        }
    }
    return norm;
}


/**
 * @brief Calculate the normalized similarity matrix
 * 
 * The function calculates the normalized similarity matrix for the given similarity matrix and diagonal degree matrix.
 * 
 * @param similarityMat - A pointer to the head of a 2D array representing the similarity matrix
 * @param diagonalDegreeMat - A pointer to the head of a 2D array representing the diagonal degree matrix
 * @param n - The number of rows and columns in the matrices
 * @return double** - A pointer to the head of a 2D array representing the normalized similarity matrix
 */
double** normalizedSimilarityMatrix(double **similarityMat, double **diagonalDegreeMat, int n)
{
    double **first;
    double **sec;
    double **diagInvSqrtMat;
    diagInvSqrtMat = diagInvSqrtMatrix(diagonalDegreeMat, n);
    first = matrixMultiply(diagInvSqrtMat, similarityMat, n, n);
    sec = matrixMultiply(first, diagInvSqrtMat, n, n);
    freeMatrix(diagInvSqrtMat, n);
    freeMatrix(first, n);
    return sec;
}

int main(int argc, char *argv[]) 
{
    char *goal = argv[1];
    FILE *file;
    char sym[] = "sym";
    char ddg[] = "ddg";
    char norm[] = "norm";

   struct vector *points;
   int numOfPoints;

   file = fopen(argv[2], "r");
   points = loadPoints(file);
   fclose(file);
   numOfPoints = countVectors(points);

    double **similarityMat = similarityMatrix(points, numOfPoints);
    if (compareStrings(goal, sym) == 1)
    {
        printMatrix(similarityMat, numOfPoints);
        
    }  else {
        double **diagonalDegreeMat = diagonalDegreeMatrix(similarityMat, numOfPoints);
        if (compareStrings(goal, ddg) == 1)
        {
            printMatrix(diagonalDegreeMat, numOfPoints);
            freeMatrix(diagonalDegreeMat, numOfPoints);
        } else {
            double **normalizedSimilarityMat = normalizedSimilarityMatrix(similarityMat, diagonalDegreeMat, numOfPoints);
            if (compareStrings(goal, norm) == 1)
            {
                printMatrix(normalizedSimilarityMat, numOfPoints);
                freeMatrix(diagonalDegreeMat, numOfPoints);
                freeMatrix(normalizedSimilarityMat, numOfPoints);
            }
        }
    }
    freeMatrix(similarityMat, numOfPoints);
    freeVector(points);
    return 0;
}
