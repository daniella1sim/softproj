#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct cord
{
    double value;
    struct cord *next;
};
struct vector
{
    struct vector *next;
    struct cord *cords;
};

struct cord* createNewCord()
{
    struct cord *headCord;
    headCord = (struct cord*)calloc(1, sizeof(struct cord));
    headCord->next = NULL;
    return headCord;
}

struct vector* createNewVector()
{
    struct vector *headVector;
    headVector = (struct vector*)calloc(1, sizeof(struct vector));
    headVector->next = NULL;
    return headVector;
}

struct cord* addNewCord(struct cord *currCord)
{
    currCord->next = (struct cord*)calloc(1, sizeof(struct cord));
    currCord = currCord->next;
    currCord->next = NULL;
    return currCord;
}

struct vector* addVector(struct vector *currVec)
{
    currVec->next = (struct vector*)calloc(1, sizeof(struct vector));
    currVec = currVec->next;
    currVec->next = NULL;
    return currVec;
}

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

struct vector* loadPoints(FILE *file)
{
    char ch ,prevCh;
    double currVal;

    struct vector *currVec, *headVec;
    struct cord *headCord, *currCord;

    headCord = (struct cord*)calloc(1, sizeof(struct cord));
    currCord = headCord;
    currCord->next = NULL;

    headVec = (struct vector*)calloc(1, sizeof(struct vector));
    currVec = headVec;
    currVec->next = NULL;

    prevCh = '\0';
    while (fscanf(file, "%lf%c", &currVal, &ch) == 2)
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

void freeMatrix(double **matrix, int n)
{
    for(int i = 0; i < n; i ++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

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

double columnSum(double** matrix, int col, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += matrix[i][col];
    }
    return sum;
}

double** matrixMultiply(double** matrix1, double** matrix2, int n)
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
            for (int k = 0; k < n; k++)
            {
                ij += matrix1[i][k] * matrix2[k][j];
            }
            matrix[i][j] = ij;
        }
    }
    return matrix;
}

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

double** similarityMatrix(struct vector *points, int numOfPoints)
{
    struct vector *currPointRow, *currPointColumn;
    double **matrix;

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

double** diagonalDegreeMatrix(double **similarityMatrix, int n)
{
    double **ddg;
    ddg = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        ddg[i] = (double*)calloc(n, sizeof(double));
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                ddg[i][j] = columnSum(similarityMatrix, j, n);
            } else {
                ddg[i][j] = 0;
            }
        }
    }
    return ddg;
}

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

double** normalizedSimilarityMatrix(double **similarityMat, double **diagonalDegreeMat, int n)
{
    double **first, **sec, **diagInvSqrtMat;
    diagInvSqrtMat = diagInvSqrtMatrix(diagonalDegreeMat, n);
    first = matrixMultiply(diagInvSqrtMat, similarityMat, n);
    sec = matrixMultiply(first, diagInvSqrtMat, n);
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
        freeMatrix(similarityMat, numOfPoints);
    }  else {
        double **diagonalDegreeMat = diagonalDegreeMatrix(similarityMat, numOfPoints);
        if (compareStrings(goal, ddg) == 1)
        {
            printMatrix(diagonalDegreeMat, numOfPoints);
            freeMatrix(similarityMat, numOfPoints);
            freeMatrix(diagonalDegreeMat, numOfPoints);
        } else {
            double **normalizedSimilarityMat = normalizedSimilarityMatrix(similarityMat, diagonalDegreeMat, numOfPoints);
            if (compareStrings(goal, norm) == 1)
            {
                printMatrix(normalizedSimilarityMat, numOfPoints);
                freeMatrix(similarityMat, numOfPoints);
                freeMatrix(diagonalDegreeMat, numOfPoints);
                freeMatrix(normalizedSimilarityMat, numOfPoints);
            }
        }
    }
    freeVector(points);
    return 1;
}